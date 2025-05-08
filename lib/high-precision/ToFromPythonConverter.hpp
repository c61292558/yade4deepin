/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef REAL_TO_FROM_PYTHON_CONVERTER_HPP
#define REAL_TO_FROM_PYTHON_CONVERTER_HPP

#include <boost/python.hpp>

#include <lib/high-precision/RealIO.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <functional>

// This is the same macro for looping over RealHP<N> levels like in lib/high-precision/RealHPConfig.hpp, but this one loops over levels exported to python (faster compilation)
#define YADE_REGISTER_HP_PYTHON_LEVELS(name) BOOST_PP_SEQ_FOR_EACH(YADE_HP_PARSE_ONE, name, YADE_MINIEIGEN_HP) // it just creates: name(1) name(2) name(3) ....

namespace forCtags {
struct ToFromPythonConverter {
}; // for ctags
}

/*************************************************************************/
/*************************       mpmath         **************************/
/*************************************************************************/

// The note at the end of http://mpmath.org/doc/current/basics.html#temporarily-changing-the-precision
// indicates that having different mpmath variables with different precision is poorly supported.
// So python conversions of RealHP<N> for different precisions is questionable.
// Not a big problem, because N>=2 is supposed to be used only in critical C++ sections where better calculations are necessary.
// And not much with python, only for debugging when necessary.
// Also, we now have directly Real and Complex types in YADE, not mpmath.mpf and mpmath.mpc anymore.
// So the prepareMpmath below is only used by RealVisitor and ComplexVisitor to export to mpmath.
template <typename Rr> struct prepareMpmath {
	static inline ::boost::python::object work()
	{
		// http://mpmath.org/doc/current/technical.html
		::boost::python::object mpmath = ::boost::python::import("mpmath"); // this code is never compiled if python3-mpmath package is unavailable.
		mpmath.attr("mp").attr("dps")  = int(std::numeric_limits<Rr>::digits10 + ::yade::math::RealHPConfig::extraStringDigits10);
		return mpmath;
	}
};

/*************************************************************************/
/*************************        Real          **************************/
/*************************************************************************/

// https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/faq/how_can_i_automatically_convert_.html
template <typename ArbitraryReal> struct ArbitraryReal_from_python {
	ArbitraryReal_from_python() { boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<ArbitraryReal>()); }
	static void* convertible(PyObject* obj_ptr)
	{
		// read from python int
		if (PyLong_CheckExact(obj_ptr) == 1) { return obj_ptr; }
		// read from python float
		if (PyFloat_CheckExact(obj_ptr) == 1) { return obj_ptr; }
		// read from all other precision RealHP<levelHP> levels.
		if (PyObject_HasAttrString(obj_ptr, "levelRealHPMethod") == 1) { return obj_ptr; }
		// read from mpmath.mpf
		if (PyObject_HasAttrString(obj_ptr, "_mpf_") == 1) { return obj_ptr; }
		// The quick way didn't work. So check if that is a string with a valid number inside.
		// This is a little more expensive. But it is used very rarely - only when user writes a python line like val="0.123123123123123123123333312312333333123123123"
		// Otherwise only Real type objects are passed around inside python scripts which does not reach this line.
		// So this approach is slow, because it's parsing a string (twice: (1) here and (2) in construct), but it's extremely rarely used.
		std::istringstream ss { ::boost::python::call_method<std::string>(obj_ptr, "__str__") };
		ArbitraryReal      r;
		ss >> r;
		// Must reach end of string .eof(), otherwise it means there were illegal characters
		return ((not ss.fail()) and (ss.eof())) ? obj_ptr : nullptr;
	}
	// NOTE: The ::boost::enable_if_c is here to make sure that ArbitraryReal_from_python::construct DOES NOT EXIST for type 'double'
	//       and produces compilation error if someone tries to use it. This is because python handles 'double' type natively.
	static typename ::boost::enable_if_c<not std::is_same<double, ArbitraryReal>::value, void>::type
	construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<ArbitraryReal>*)(data))->storage.bytes;
		new (storage) ArbitraryReal;
		ArbitraryReal* val     = (ArbitraryReal*)storage;
		bool           success = false;
		// read from python int
		if (PyLong_CheckExact(obj_ptr) == 1) {
			*val    = boost::python::extract<long int>(obj_ptr)();
			success = true;
		}
		// read from python float
		if ((not success) and (PyFloat_CheckExact(obj_ptr) == 1)) {
			*val    = boost::python::extract<double>(obj_ptr)();
			success = true;
		}
		// RealVisitor : handle all other high-precision types RealHP<levelHP>
		if ((not success) and (PyObject_HasAttrString(obj_ptr, "levelRealHPMethod") == 1)) {
			int objLevelHP = ::boost::python::call_method<int>(obj_ptr, "levelRealHPMethod");
			switch (objLevelHP) {
#define CASE_LEVEL_HP(levelHP)                                                                                                                                 \
	case levelHP: {                                                                                                                                        \
		if (objLevelHP == ::yade::math::levelOfRealHP<ArbitraryReal>) {                                                                                \
			throw std::runtime_error(__FILE__ " :  ArbitraryReal_from_python::construct objLevelHP == levelOfRealHP<ArbitraryReal> error.");       \
		} else {                                                                                                                                       \
			*val = static_cast<ArbitraryReal>(boost::python::extract<::yade::math::RealHP<levelHP>>(obj_ptr)());                                   \
		}                                                                                                                                              \
	} break;
				YADE_REGISTER_HP_PYTHON_LEVELS(CASE_LEVEL_HP)
#undef CASE_LEVEL_HP
				default: throw std::runtime_error(__FILE__ " :  ArbitraryReal_from_python::construct switch default case objLevelHP error.");
			}
			success = true;
		}
		// No need to check this: ↘ because taking from mpmath.mpf is the same as taking from a string.
		if ((not success) /* or (PyObject_HasAttrString(obj_ptr, "_mpf_") == 1) */) {
			prepareMpmath<ArbitraryReal>::work();
			std::istringstream ss { ::boost::python::call_method<std::string>(obj_ptr, "__str__") };
			*val    = ::yade::math::fromStringRealHP<ArbitraryReal>(ss.str());
			success = true;
		}
		/* Not needed in fact. Last string conversion covered all possible cases.
		if (not success) {
			throw std::runtime_error(__FILE__ " :  ArbitraryReal_from_python::construct failed to construct a Real type from given type.");
		}
		*/
		data->convertible = storage;
	}
};

/*************************************************************************/
/*************************       Complex        **************************/
/*************************************************************************/

// TODO: the convertible(…) and construct(…) for Complex type below are in fact extremely slow, because they work on strings.
//       If we ever need fast high precision Complex numbers in YADE, then these two functions should be rewritten
//       to work in similar manner as ArbitraryReal_from_python::convertible(…) and construct(…)
// https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/faq/how_can_i_automatically_convert_.html
template <typename ArbitraryComplex> struct ArbitraryComplex_from_python {
	ArbitraryComplex_from_python() { boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<ArbitraryComplex>()); }
	static void* convertible(PyObject* obj_ptr)
	{
		// only python complex or mpmath.mpc(…) objects are supoprted. Strings are not parsed.
		// However a simple workaround is to write mpmath.mpc("1.211213123123123123123123123","-124234234.111")
		PyComplex_AsCComplex(obj_ptr);
		if (PyErr_Occurred() == nullptr) return obj_ptr;
		PyErr_Clear();
		return nullptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		prepareMpmath<typename ArbitraryComplex::value_type>::work();
		std::istringstream ss_real { ::boost::python::call_method<std::string>(
			::boost::python::expect_non_null(PyObject_GetAttrString(obj_ptr, "real")), "__str__") };
		std::istringstream ss_imag { ::boost::python::call_method<std::string>(
			::boost::python::expect_non_null(PyObject_GetAttrString(obj_ptr, "imag")), "__str__") };
		void*              storage = ((boost::python::converter::rvalue_from_python_storage<ArbitraryComplex>*)(data))->storage.bytes;
		new (storage) ArbitraryComplex;
		ArbitraryComplex*                     val = (ArbitraryComplex*)storage;
		typename ArbitraryComplex::value_type re { 0 }, im { 0 };
		// ensure that "nan" "inf" are read correctly
		re                = ::yade::math::fromStringRealHP<typename ArbitraryComplex::value_type>(ss_real.str());
		im                = ::yade::math::fromStringRealHP<typename ArbitraryComplex::value_type>(ss_imag.str());
		*val              = ArbitraryComplex(re, im); // must explicitly call the constructor, static_cast won't work.
		data->convertible = storage;
	}
};

/*************************************************************************/
/*************************   RealHP + python    **************************/
/*************************************************************************/

namespace yade {
namespace math {
	namespace detail {
		template <int, template <int, bool> class>
		class ScopeHP { // separate class is needed to act as python scope identifier. Might become useful later.
		};

		template <int N, template <int, bool> class RegisterHPClass> void registerInScope(bool createInternalScopeHP)
		{
			boost::python::scope topScope;
			if (createInternalScopeHP) {
				// This creates internal python scope HP1 or HP2 or HP3 and so on. In each of them are the same math functions with respective precisions.
				// The main point is that all math functions, and minieigen classes are accessible e.g. for RealHP<4>, via:
				//    yade.math.HP4.sin(10)
				//    yade.minieigenHP.HP4.Vector3r(1,2,3)
				// The original RealHP<1> are present in two places:
				//    yade.math.sin(10)                             yade.math.HP1.sin(10)
				//    yade.minieigenHP.Vector3r(1,2,3)              yade.minieigenHP.HP1.Vector3r(1,2,3)
				std::string name = "HP" + boost::lexical_cast<std::string>(N); // scope name: "HP1", "HP2", etc
				if (PyObject_HasAttrString(topScope.ptr(), name.c_str()) == true) {
					// If the scope already exists, then use it to add more objects to it.
					RegisterHPClass<N, true>::work(topScope, boost::python::scope(topScope.attr(name.c_str())));
				} else {
					boost::python::scope HPn = boost::python::class_<ScopeHP<N, RegisterHPClass>>(name.c_str());
					RegisterHPClass<N, true>::work(topScope, HPn);
				}
			} else {
				// this one puts the 'base precision' RealHP<1> math functions in the top scope of this python module. They are duplicated inside HP1.
				// Not sure which place is more convenient to use. Maybe both.
				RegisterHPClass<N, false>::work(topScope, topScope);
			}
		}

		template <template <int, bool> class RegisterHPClass> struct WorkaroundClangCompiler {
			template <typename N1> void operator()(N1) { registerInScope<N1::value, RegisterHPClass>(true); }
		};

		// this loop registers python functions from Range by calling the provided RegisterHPClass<int,bool>::work( , ); inside registerInScope above.
		template <typename Range, template <int, bool> class RegisterHPClass> void registerLoopForHPn()
		{
			registerInScope<1, RegisterHPClass>(false);
			boost::mpl::for_each<Range>(WorkaroundClangCompiler<RegisterHPClass>());
		}
	}
}
}

/*************************************************************************/
/************************* minieigenHP → string **************************/
/*************************************************************************/

// these are used by py/high-precision/minieigen/visitors.hpp
namespace yade {
namespace minieigenHP {
	template <typename Rr, typename boost::disable_if_c<::yade::math::isComplex<Rr>, int>::type = 0, int Level = ::yade::math::levelOfRealHP<Rr>>
	inline std::string numToStringHP(const Rr& num)
	{
		using R = typename std::decay<Rr>::type;
		static_assert(std::is_same<R, ::yade::RealHP<Level>>::value, "RealHP problem here, please file a bug report.");
		constexpr bool isPythonPrecisionEnough = std::is_same<double, R>::value or std::is_same<float, R>::value;
		if (isPythonPrecisionEnough) {
			return ::yade::math::toStringHP(num);
		} else {
			// The only way to make sure that it is copy-pasteable to/from python without loss of precision is to put the numbers inside "…"
			return "\"" + ::yade::math::toStringHP(num) + "\"";
		}
	}

	template <typename Cc, typename boost::enable_if_c<::yade::math::isComplex<Cc>, int>::type = 0, int Level = ::yade::math::levelOfComplexHP<Cc>>
	inline std::string numToStringHP(const Cc& num)
	{
		using C = typename std::decay<Cc>::type;
		using R = typename C::value_type;
		static_assert(std::is_same<C, ::yade::ComplexHP<Level>>::value, "ComplexHP problem here, please file a bug report.");
		constexpr bool isPythonPrecisionEnough = std::is_same<double, R>::value or std::is_same<float, R>::value;
		std::string    ret;
		if (num.real() != 0 && num.imag() != 0) {
			if (isPythonPrecisionEnough) {
				// don't add "+" in the middle if imag is negative and will start with "-"
				return numToStringHP(num.real()) + (num.imag() > 0 ? "+" : "") + numToStringHP(num.imag()) + "j";
			} else {
				// make sure it is copy-pasteable without loss of precision
				return "Complex(" + numToStringHP(num.real()) + "," + numToStringHP(num.imag()) + ")";
			}
		}
		// only imaginary is non-zero: skip the real part
		if (num.imag() != 0) {
			if (isPythonPrecisionEnough) {
				return numToStringHP(num.imag()) + "j";
			} else {
				return "Complex(\"0\"," + numToStringHP(num.imag()) + ")";
			}
		}
		if (isPythonPrecisionEnough) {
			return numToStringHP(num.real());
		} else {
			return "Complex(" + numToStringHP(num.real()) + ",\"0\")";
		}
	}

	inline std::string numToStringHP(const int& num) { return ::boost::lexical_cast<::std::string>(num); } // ignore padding for now.

} // namespace minieigenHP
} // namespace yade

#endif
