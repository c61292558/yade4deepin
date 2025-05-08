/*************************************************************************
*  2025      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;

#include <py/high-precision/minieigen/RealVisitor.hpp>
template <int N, bool> struct Helper {
	static void expose_math_Real(bool /* notDuplicate */, const py::scope& /* topScope */) {};
};

template <int N> struct Helper<N, true> {
	static void expose_math_Real(bool notDuplicate, const py::scope& topScope)
	{
		if (notDuplicate) {
			py::class_<RealHP<N>>(
			        "Real",
			        "The Real type.", // TODO: add more description
			        py::init<>())
			        .def(RealVisitor<RealHP<N>>());
		} else {
			py::scope().attr("Real") = topScope.attr("Real");
		}
	};
};

template <int N> void expose_math_Real(bool notDuplicate, const py::scope& topScope)
{
	// We have almost switched to C++17 which supports 'if constexpr'. But some compile targets:
	//   make_18_04, make_18_04_nogui, make_suse15, mk_bll_arm64, mk_bll_ppc64el, mk_bookw_arm64 use C++14
	// still use C++14, so I cannot use 'if constexpr' and instead have to use a clunky template
	Helper<N, (std::numeric_limits<RealHP<N>>::digits10 >= 18)>::expose_math_Real(notDuplicate, topScope);
}

// explicit instantination - tell compiler to produce a compiled version of expose_converters (it is faster when done in parallel in .cpp files)
YADE_HP_PYTHON_REGISTER(expose_math_Real)

// Call TEMPLATE_CREATE_LOGGER for all high-precision instances.
// To see debug messages from constructors use command:
//    yade.log.setLevel("RealVisitor.hpp", 4)
// But uncomment all calls to LOG_* in RealVisitor.hpp first.
#define YADE_HP_PY_REAL_LOGGER(r, name, levelHP) TEMPLATE_CREATE_LOGGER(name<RealHP<levelHP>>);
#define YADE_HP_PYTHON_REAL_TEMPLATE_CREATE_LOGGER(name)                                                                                                       \
	BOOST_PP_SEQ_FOR_EACH(YADE_HP_PY_REAL_LOGGER, name, YADE_MINIEIGEN_HP) // instatinate templates for name<1>, name<2>, etc …
YADE_HP_PYTHON_REAL_TEMPLATE_CREATE_LOGGER(RealVisitor);
