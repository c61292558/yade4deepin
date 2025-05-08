/*************************************************************************
*  2025      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#include "common.hpp"
#include <lib/high-precision/MathComplexFunctions.hpp>

#include <lib/base/LoggingUtils.hpp>

template <typename ComplexT, int Level = ::yade::math::levelOfComplexHP<ComplexT>> class ComplexVisitor : public py::def_visitor<ComplexVisitor<ComplexT>> {
	typedef typename ComplexT::value_type RealT;

public:
	template <class PyClass> void visit(PyClass& cl) const
	{
		// one argument
		cl.def("__init__", py::make_constructor(&ComplexVisitor::fromObj, py::default_call_policies(), (py::arg("obj")))) // mpmath support
#if (YADE_REAL_BIT >= 80)
		        .def("__init__", py::make_constructor(&ComplexVisitor::fromComplexT, py::default_call_policies(), (py::arg("z"))))
		        .def("__init__", py::make_constructor(&ComplexVisitor::fromRealT, py::default_call_policies(), (py::arg("z"))))
#endif
		        .def("__init__", py::make_constructor(&ComplexVisitor::fromComplex, py::default_call_policies(), (py::arg("z"))))
		        .def("__init__", py::make_constructor(&ComplexVisitor::fromDouble, py::default_call_policies(), (py::arg("d"))))
		        .def("__init__", py::make_constructor(&ComplexVisitor::fromInt, py::default_call_policies(), (py::arg("i"))))
		        .def("__init__", py::make_constructor(&ComplexVisitor::fromString, py::default_call_policies(), (py::arg("str"))))
		        // two arguments, re, im.
		        .def("__init__",
		             py::make_constructor(&ComplexVisitor::from2Objs, py::default_call_policies(), (py::arg("re"), py::arg("im")))) // mpmath support
#if (YADE_REAL_BIT >= 80)
		        .def("__init__",
		             py::make_constructor(&ComplexVisitor::from2RealT, py::default_call_policies(), (py::arg("re"), py::arg("im")))) // mpmath support
#endif
		        .def("__init__", py::make_constructor(&ComplexVisitor::from2Doubles, py::default_call_policies(), (py::arg("a"), py::arg("b"))))
		        .def("__init__", py::make_constructor(&ComplexVisitor::from2Ints, py::default_call_policies(), (py::arg("i"), py::arg("j"))))
		        .def("__init__", py::make_constructor(&ComplexVisitor::from2Strings, py::default_call_policies(), (py::arg("str1"), py::arg("str2"))))
		        .def_pickle(ComplexPickle())
		        // operators
		        .def(py::self + py::self)
		        .def(py::self + py::other<RealT>())
		        .def(py::other<RealT>() + py::self)
		        .def(py::self += py::self)
		        .def(py::self += py::other<RealT>())
		        .def(py::self - py::self)
		        .def(py::self - py::other<RealT>())
		        .def(py::other<RealT>() - py::self)
		        .def(py::self -= py::self)
		        .def(py::self -= py::other<RealT>())
		        .def(py::self * py::self)
		        .def(py::self * py::other<RealT>())
		        .def(py::other<RealT>() * py::self)
		        .def(py::self *= py::self)
		        .def(py::self *= py::other<RealT>())
		        .def(py::self / py::self)
		        .def(py::self / py::other<RealT>())
		        .def(py::other<RealT>() / py::self)
		        .def(py::self /= py::self)
		        .def(py::self /= py::other<RealT>())
		        .def("__eq__", &ComplexVisitor::__eq__)
		        .def("__ne__", &ComplexVisitor::__ne__)
		        // specials
		        .add_property("real", &ComplexVisitor::real)
		        .add_property("imag", &ComplexVisitor::imag)
		        .def("levelComplexHPMethod", &ComplexVisitor::levelHP)
		        .add_property("levelHP", &ComplexVisitor::levelHP)
		        .def("__pos__", &ComplexVisitor::__pos__)
		        .def("__neg__", &ComplexVisitor::__neg__)
		        .def("__abs__", &ComplexVisitor::__abs__)
		        .def("__complex__", &ComplexVisitor::toComplex)
		        .def("__str__", &ComplexVisitor::__str__)
		        .def("__repr__", &ComplexVisitor::__repr__)
		        .add_property("_mpc_", &ComplexVisitor::_mpc_);
	}

private:
	static int       levelHP(const ComplexT&) { return Level; }
	static ComplexT* fromObj(const py::object& obj) // mpmath support
	{
		RealT     re  = ::yade::math::fromStringRealHP<RealT>(py::extract<std::string>(obj.attr("real").attr("__str__")()));
		RealT     im  = ::yade::math::fromStringRealHP<RealT>(py::extract<std::string>(obj.attr("imag").attr("__str__")()));
		ComplexT* ret = new ComplexT(re, im);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from slow mpmath.mpc(…). Better to use " << ::yade::minieigenHP::numToStringHP(*ret) << " instead.");
		return ret;
	}
	static ComplexT* fromString(const std::string& str)
	{
		ComplexT* ret = new ComplexT(::yade::math::fromStringRealHP<RealT, Level>(str));
		// Parsing in fromStringComplexHP is a bit strange - requires argument like "(1,2)" which is the same as from2Strings.
		// So lets use the above from Real number instead.
		// Can be improved later.
		//ComplexT* ret = new ComplexT(::yade::math::fromStringComplexHP<ComplexT, Level>(str));
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from string, got: " << ::yade::minieigenHP::numToStringHP(*ret) << ", note: this function needs a fix if is to be used more often.");
		return ret;
	}
	static ComplexT* fromComplexT(const ComplexT& d)
	{
		ComplexT* ret = new ComplexT(d);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from ComplexT, got: " << ::yade::minieigenHP::numToStringHP(*ret) << ", good.");
		return ret;
	}
	static ComplexT* fromRealT(const RealT& d)
	{
		ComplexT* ret = new ComplexT(d);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from RealT, got: " << ::yade::minieigenHP::numToStringHP(*ret) << ", good.");
		return ret;
	}
	static ComplexT* fromDouble(const double& d)
	{
		ComplexT* ret = new ComplexT(d);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from python float, got: " << ::yade::minieigenHP::numToStringHP(*ret));
		return ret;
	}
	static ComplexT* fromComplex(const ::std::complex<double>& d)
	{
		ComplexT* ret = new ComplexT(d);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from python complex, got: " << ::yade::minieigenHP::numToStringHP(*ret));
		return ret;
	}
	static ComplexT* fromInt(const long int& i)
	{
		ComplexT* ret = new ComplexT(i);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from python int, got: " << ::yade::minieigenHP::numToStringHP(*ret));
		return ret;
	}
	static ComplexT* from2Objs(const py::object& re_obj, const py::object& im_obj) // mpmath support
	{
		RealT     re  = ::yade::math::fromStringRealHP<RealT>(py::extract<std::string>(re_obj.attr("__str__")()));
		RealT     im  = ::yade::math::fromStringRealHP<RealT>(py::extract<std::string>(im_obj.attr("__str__")()));
		ComplexT* ret = new ComplexT(re, im);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from two slow mpmath arguments. Better to use " << ::yade::minieigenHP::numToStringHP(*ret) << " instead.");
		return ret;
	}
	static ComplexT* from2Strings(const std::string& str1, const std::string& str2)
	{
		ComplexT* ret = new ComplexT(::yade::math::fromStringRealHP<RealT, Level>(str1), ::yade::math::fromStringRealHP<RealT, Level>(str2));
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from two strings, got: " << ::yade::minieigenHP::numToStringHP(*ret));
		return ret;
	}
	static ComplexT* from2RealT(const RealT& a, const RealT& b)
	{
		ComplexT* ret = new ComplexT(a, b);
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from two python RealT, got: " << ::yade::minieigenHP::numToStringHP(*ret));
		return ret;
	}
	static ComplexT* from2Doubles(const double& a, const double& b)
	{
		ComplexT* ret = new ComplexT(RealT(a), RealT(b));
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from two python floats, got: " << ::yade::minieigenHP::numToStringHP(*ret));
		return ret;
	}
	static ComplexT* from2Ints(const long int& i, const long int& j)
	{
		ComplexT* ret = new ComplexT(RealT(i), RealT(j));
		//LOG_TIMED_4_INFO(1s, "Constructing Complex from two long int arguments, got: " << ::yade::minieigenHP::numToStringHP(*ret));
		return ret;
	}

	struct ComplexPickle : py::pickle_suite {
		static py::tuple getinitargs(const ComplexT& x) { return py::make_tuple(x.real(), x.imag()); }
	};
	static bool                 __eq__(const ComplexT& u, const ComplexT& v) { return u == v; }
	static bool                 __ne__(const ComplexT& u, const ComplexT& v) { return !__eq__(u, v); }
	static RealT                real(const ComplexT& a) { return a.real(); }
	static RealT                imag(const ComplexT& a) { return a.imag(); }
	static ComplexT             __pos__(const ComplexT& a) { return +a; }
	static ComplexT             __neg__(const ComplexT& a) { return -a; }
	static RealT                __abs__(const ComplexT& a) { return ::yade::math::abs(a); }
	static std::complex<double> toComplex(const ComplexT& a) { return std::complex<double>(static_cast<double>(a.real()), static_cast<double>(a.imag())); }

	static std::string __str__(const py::object& obj)
	{
		const ComplexT self = py::extract<ComplexT>(obj)();
		if (self.imag() >= 0) {
			return "(" + ::yade::math::toStringHP(self.real()) + "+" + ::yade::math::toStringHP(self.imag()) + "j)";
		} else {
			return "(" + ::yade::math::toStringHP(self.real()) + "-" + ::yade::math::toStringHP(-self.imag()) + "j)";
		}
	}
	static std::string __repr__(const py::object& obj)
	{
		const ComplexT self = py::extract<ComplexT>(obj)();
		return std::string(object_class_name(obj) + "(") + ::yade::minieigenHP::numToStringHP(self.real()) + ","
		        + ::yade::minieigenHP::numToStringHP(self.imag()) + ")";
	}
	static py::tuple _mpc_(const py::object& obj)
	{
		const ComplexT val = py::extract<ComplexT>(obj)();

		//LOG_TIMED_3_WARN(2s, "Exporting Complex " << ::yade::minieigenHP::numToStringHP(val) << " to slow mpmath. Try to avoid this by working directly with Complex(…)");

		std::stringstream ss_real {};
		std::stringstream ss_imag {};
		ss_real << ::yade::math::toStringHP<RealT>(val.real());
		ss_imag << ::yade::math::toStringHP<RealT>(val.imag());
		py::object mpmath = prepareMpmath<RealT>::work();
		py::object result = mpmath.attr("mpc")(ss_real.str(), ss_imag.str());
		return py::extract<py::tuple>(result.attr("_mpc_"));
	}
	DECLARE_LOGGER;
};
