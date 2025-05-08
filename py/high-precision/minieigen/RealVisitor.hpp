/*************************************************************************
*  2025      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#include "common.hpp"

#include <lib/base/LoggingUtils.hpp>

template <typename RealT, int Level = ::yade::math::levelOfRealHP<RealT>> class RealVisitor : public py::def_visitor<RealVisitor<RealT>> {
public:
	template <class PyClass> void visit(PyClass& cl) const
	{
		cl.def("__init__", py::make_constructor(&RealVisitor::fromObj, py::default_call_policies(), (py::arg("obj")))) // mpmath support
#if (YADE_REAL_BIT >= 80)
		        .def("__init__", py::make_constructor(&RealVisitor::fromRealT, py::default_call_policies(), (py::arg("z"))))
#endif
		        .def("__init__", py::make_constructor(&RealVisitor::fromDouble, py::default_call_policies(), (py::arg("d"))))
		        .def("__init__", py::make_constructor(&RealVisitor::fromInt, py::default_call_policies(), (py::arg("i"))))
		        .def("__init__", py::make_constructor(&RealVisitor::fromString, py::default_call_policies(), (py::arg("str"))))
		        .def_pickle(RealPickle())
		        // operators
		        .def(py::self + py::self)
		        .def(py::self + py::other<long int>())
		        .def(py::self + py::other<double>())
		        .def(py::other<long int>() + py::self)
		        .def(py::other<double>() + py::self)
		        .def(py::self += py::self)
		        .def(py::self += py::other<long int>())
		        .def(py::self += py::other<double>())
		        .def(py::self - py::self)
		        .def(py::self - py::other<long int>())
		        .def(py::self - py::other<double>())
		        .def(py::other<long int>() - py::self)
		        .def(py::other<double>() - py::self)
		        .def(py::self -= py::self)
		        .def(py::self -= py::other<long int>())
		        .def(py::self -= py::other<double>())
		        .def(py::self * py::self)
		        .def(py::self * py::other<long int>())
		        .def(py::self * py::other<double>())
		        .def(py::other<long int>() * py::self)
		        .def(py::other<double>() * py::self)
		        .def(py::self *= py::self)
		        .def(py::self *= py::other<long int>())
		        .def(py::self *= py::other<double>())
		        .def(py::self / py::self)
		        .def(py::self / py::other<long int>())
		        .def(py::self / py::other<double>())
		        .def(py::other<long int>() / py::self)
		        .def(py::other<double>() / py::self)
		        .def(py::self /= py::self)
		        .def(py::self /= py::other<long int>())
		        .def(py::self /= py::other<double>())
		        .def(py::self < py::self)
		        .def(py::self < py::other<long int>())
		        .def(py::self < py::other<double>())
		        .def(py::other<long int>() < py::self)
		        .def(py::other<double>() < py::self)
		        .def(py::self > py::self)
		        .def(py::self > py::other<long int>())
		        .def(py::self > py::other<double>())
		        .def(py::other<long int>() > py::self)
		        .def(py::other<double>() > py::self)
		        .def(py::self <= py::self)
		        .def(py::self <= py::other<long int>())
		        .def(py::self <= py::other<double>())
		        .def(py::other<long int>() <= py::self)
		        .def(py::other<double>() <= py::self)
		        .def(py::self >= py::self)
		        .def(py::self >= py::other<long int>())
		        .def(py::self >= py::other<double>())
		        .def(py::other<long int>() >= py::self)
		        .def(py::other<double>() >= py::self)
		        .def(py::self == py::other<long int>())
		        .def(py::self == py::other<double>())
		        .def(py::other<long int>() == py::self)
		        .def(py::other<double>() == py::self)
		        .def(py::self != py::other<long int>())
		        .def(py::self != py::other<double>())
		        .def(py::other<long int>() != py::self)
		        .def(py::other<double>() != py::self)
		        .def("__eq__", &RealVisitor::__eq__)
		        .def("__ne__", &RealVisitor::__ne__)
		        // specials
		        .add_property("real", &RealVisitor::real)
		        .add_property("imag", &RealVisitor::imag)
		        .def("levelRealHPMethod", &RealVisitor::levelHP)
		        .add_property("levelHP", &RealVisitor::levelHP)
		        .def("__pow__", &RealVisitor::__powInt__)
		        .def("__pow__", &RealVisitor::__powReal__)
		        .def("sqrt", &RealVisitor::__sqrt__)
		        .def("__mod__", &RealVisitor::__mod__)
		        .def("__round__", &RealVisitor::__round__)
		        .def("__round__", &RealVisitor::__roundNdigits__)
		        .def("__pos__", &RealVisitor::__pos__)
		        .def("__neg__", &RealVisitor::__neg__)
		        .def("__abs__", &RealVisitor::__abs__)
		        .def("__float__", &RealVisitor::toFloat)
		        .def("__int__", &RealVisitor::toInt)
		        .def("__str__", &RealVisitor::__str__)
		        .def("__repr__", &RealVisitor::__repr__)
		        .add_property("_mpf_", &RealVisitor::_mpf_);
	}

private:
	static int                levelHP(const RealT&) { return Level; }
	static inline std::string fullStr(const RealT& val) { return "Real(" + ::yade::minieigenHP::numToStringHP(val) + ")"; }
	static RealT*             fromObj(const py::object& obj) // mpmath support
	{
		RealT* ret = new RealT(::yade::math::fromStringRealHP<RealT>(py::extract<std::string>(obj.attr("real").attr("__str__")())));
		//LOG_TIMED_4_INFO(1s, "Constructing Real from slow mpmath.mpf(…). Better to use " << fullStr(*ret) << " instead.");
		return ret;
	}
	static RealT* fromString(const std::string& str)
	{
		RealT* ret = new RealT(::yade::math::fromStringRealHP<RealT, Level>(str));
		//LOG_TIMED_4_INFO(1s, "Constructing Real from string, got: " << fullStr(*ret));
		return ret;
	}
	static RealT* fromRealT(const RealT& d)
	{
		RealT* ret = new RealT(d);
		//LOG_TIMED_4_INFO(1s, "Constructing Real from RealT, got: " << fullStr(*ret) << ", good.");
		return ret;
	}
	static RealT* fromDouble(const double& d)
	{
		RealT* ret = new RealT(d);
		//LOG_TIMED_4_INFO(1s, "Constructing Real from python float, got: " << fullStr(*ret));
		return ret;
	}
	static RealT* fromInt(const long int& i)
	{
		RealT* ret = new RealT(i);
		//LOG_TIMED_4_INFO(1s, "Constructing Real from python long int, got: " << fullStr(*ret));
		return ret;
	}

	struct RealPickle : py::pickle_suite {
		static py::tuple getinitargs(const RealT& x) { return py::make_tuple(x); }
	};
	static bool  __eq__(const RealT& u, const RealT& v) { return u == v; }
	static bool  __ne__(const RealT& u, const RealT& v) { return !__eq__(u, v); }
	static RealT real(const RealT& a) { return a; }
	static RealT imag(const RealT&) { return 0; }
	static RealT __powInt__(const RealT& a, const long int& i) { return ::yade::math::pow(a, i); }
	static RealT __powReal__(const RealT& a, const RealT& b) { return ::yade::math::pow(a, b); }
	static RealT __sqrt__(const RealT& a) { return ::yade::math::sqrt(a); }
	static RealT __mod__(const RealT& a, const RealT& b) { return ::yade::math::fmod(a, b); }
	static RealT __round__(const RealT& a) { return ::yade::math::rint(a); }
	static RealT __roundNdigits__(const RealT& a, const int& p)
	{
		RealT mult = ::yade::math::pow(RealT(10), RealT(p));
		return ::yade::math::rint(a * mult) / mult;
	}
	static RealT  __pos__(const RealT& a) { return +a; }
	static RealT  __neg__(const RealT& a) { return -a; }
	static RealT  __abs__(const RealT& a) { return ::yade::math::abs(a); }
	static double toFloat(const RealT& a)
	{
		//LOG_TIMED_3_WARN(2s, "Converting " << fullStr(a) << " to python float and losing precision! Try to work with Real(…) instead.\n");
		return static_cast<double>(a);
	}
	static long int toInt(const RealT& a)
	{
		//LOG_TIMED_3_WARN(2s, "Converting " << fullStr(a) << " to python long int and losing precision! Try to work with Real(…) instead.\n");
		return static_cast<long int>(a);
	}

	static std::string __str__(const py::object& obj)
	{
		const RealT self = py::extract<RealT>(obj)();
		return ::yade::math::toStringHP(self);
	}
	static std::string __repr__(const py::object& obj)
	{
		const RealT self = py::extract<RealT>(obj)();
		return std::string(object_class_name(obj) + "(") + ::yade::minieigenHP::numToStringHP(self) + ")";
	}
	static py::tuple _mpf_(const py::object& obj)
	{
		const RealT val = py::extract<RealT>(obj)();

		//LOG_TIMED_3_WARN(2s, "Exporting Real " << fullStr(val) << " to slow mpmath. Try to avoid this by working directly with Real(…)");

		py::object mpmath = prepareMpmath<RealT>::work();
		if (yade::math::isnan(val)) { // mpmath does not tolerate '-nan' values. So any NaN will be a "nan".
			py::object result = mpmath.attr("mpf")("nan");
			return py::extract<py::tuple>(result.attr("_mpf_"));
		} else { // other numbers are normally converted
			py::object result = mpmath.attr("mpf")(::yade::math::toStringHP<RealT>(val));
			return py::extract<py::tuple>(result.attr("_mpf_"));
		}
	}
	DECLARE_LOGGER;
};
