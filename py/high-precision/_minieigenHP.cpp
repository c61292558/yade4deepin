/*************************************************************************
*  2009-2012 © Václav Šmilauer                                           *
*  2019        Janek Kozicki                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// When yade uses high-precision number as Real type the usual (old) 'import minieigen'
// has to be replaced with a minieigenHP library which uses exacly the same number of decimal places
// as yade is using everywhere else. Please note that such precision can be very arbitrary, because cpp_bin_float
// or mpfr take it as a compile-time argument. Hence such minieigenHP cannot be separately precompiled as a package.
// Though it could be precompiled for some special types such as boost::multiprecision::float128

#include <lib/base/Logging.hpp>
#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/RealHPConfig.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
#include <lib/pyutil/doc_opts.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <py/high-precision/_ExposeStorageOrdering.hpp>
#include <sstream>

using namespace ::yade::MathEigenTypes;
#include <py/high-precision/minieigen/converters.hpp>
#include <py/high-precision/minieigen/expose.hpp>
#include <py/high-precision/minieigen/visitors.hpp>

CREATE_CPP_LOCAL_LOGGER("_minieigenHP.cpp")

template <int N, bool> struct RegHelper {
	static void work(const py::scope& /* topScope */) {};
};

template <int N> struct RegHelper<N, true> {
	static void work(const py::scope& topScope)
	{
		py::scope top(topScope);
		ArbitraryComplex_from_python<ComplexHP<N>>();
		ArbitraryReal_from_python<RealHP<N>>();
	};
};

template <int N, bool registerConverters> struct RegisterEigenHP {
	// registerConverters is because C++ ↔ python converers need to be registered. In one of them (parent or child) it has to be skipped to avoid duplicate registrations.
	// the classes (e.g. HP8.Vector3r ↔ Vector3rHP<8>) have to be exposed always.
	static void work(const py::scope& topScope, const py::scope& scopeHP)
	{
		constexpr bool notDuplicate = not((N == 1) and registerConverters);
		RegHelper<N, (/* notDuplicate and */ (std::numeric_limits<RealHP<N>>::digits10 >= 18))>::work(topScope);
		py::scope HPn(scopeHP);
		if (std::numeric_limits<RealHP<N>>::digits10 >= 18) {
			expose_math_Complex<N>(notDuplicate, topScope);
			expose_math_Real<N>(notDuplicate, topScope);
		} else { // Real is python float type
			py::scope().attr("Real")    = py::import("__main__").attr("__dict__")["__builtins__"].attr("float");
			py::scope().attr("Complex") = py::import("__main__").attr("__dict__")["__builtins__"].attr("complex");
		}

		expose_converters<N>(notDuplicate, topScope); // in _ExposeConverters.cpp

#ifdef EIGEN_DONT_ALIGN
		py::scope().attr("vectorize") = false;
#else
		py::scope().attr("vectorize") = true;
		// when we use vectorization the Vector3r is AlignedVector3, so we need to register converter from plain old Vector3<Real> so that other functions can accept it as an argument
		if (notDuplicate) {
			custom_VectorAnyAny_from_sequence<Eigen::Matrix<RealHP<N>, 3, 1>>();
			py::class_<Eigen::Matrix<RealHP<N>, 3, 1>>(
			        "Vector3na",
			        "3-dimensional non-aligned float vector; same as :obj:`Vector3`, but with alignment (``Eigen::AlignedVector3``).\n\nSupported "
			        "operations (``f`` if a float/int, ``v`` is a Vector3): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, "
			        "``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``, plus operations with ``Matrix3`` and ``Quaternion``.\n\nImplicit conversion "
			        "from sequence (list, tuple, ...) of 3 floats.\n\nStatic attributes: ``Zero``, ``Ones``, ``UnitX``, ``UnitY``, ``UnitZ``.",
			        py::init<>())
			        .def(VectorVisitor<Eigen::Matrix<RealHP<N>, 3, 1>>());
		} else {
			py::scope().attr("Vector3na") = topScope.attr("Vector3na");
		}
#endif

		expose_vectors1<N>(notDuplicate, topScope);
		expose_vectors2<N>(notDuplicate, topScope);
		expose_matrices1<N>(notDuplicate, topScope); // must come after vectors
		expose_matrices2<N>(notDuplicate, topScope); // must come after vectors
		expose_complex1<N>(notDuplicate, topScope);
		expose_complex2<N>(notDuplicate, topScope);
		expose_quaternion<N>(notDuplicate, topScope);
		expose_boxes<N>(notDuplicate, topScope);
	}
};

BOOST_PYTHON_MODULE(_minieigenHP)
try {
	YADE_SET_DOCSTRING_OPTS;

	py::scope().attr("__doc__") = "miniEigen is wrapper for a small part of the `Eigen <http://eigen.tuxfamily.org>`_ library. Refer to its documentation "
	                              "for details. All classes in this module support pickling.";

	::yade::math::detail::registerLoopForHPn<::yade::math::RealHPConfig::SupportedByMinieigen, RegisterEigenHP>();

} catch (...) {
	LOG_FATAL("Importing this module caused an exception and this module is in an inconsistent state now.");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__);
	boost::python::handle_exception();
	throw;
}
