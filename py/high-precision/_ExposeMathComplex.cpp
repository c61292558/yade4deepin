/*************************************************************************
*  2025      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;

#include <py/high-precision/minieigen/ComplexVisitor.hpp>
template <int N> void expose_math_Complex(bool notDuplicate, const py::scope& topScope)
{
	if (std::numeric_limits<RealHP<N>>::digits10 >= 18) {
		if (notDuplicate) {
			py::class_<ComplexHP<N>>(
			        "Complex",
			        "The Complex type.", // TODO: add more description
			        py::init<>())
			        .def(ComplexVisitor<ComplexHP<N>>());
		} else {
			py::scope().attr("Complex") = topScope.attr("Complex");
		}
	}
}

// explicit instantination - tell compiler to produce a compiled version of expose_converters (it is faster when done in parallel in .cpp files)
YADE_HP_PYTHON_REGISTER(expose_math_Complex)

// Call TEMPLATE_CREATE_LOGGER for all high-precision instances.
// To see debug messages from constructors use command:
//    yade.log.setLevel("ComplexVisitor.hpp", 4)
// But uncomment all calls to LOG_* in ComplexVisitor.hpp first.
#define YADE_HP_PY_COMPLEX_LOGGER(r, name, levelHP) TEMPLATE_CREATE_LOGGER(name<ComplexHP<levelHP>>);
#define YADE_HP_PYTHON_COMPLEX_TEMPLATE_CREATE_LOGGER(name)                                                                                                    \
	BOOST_PP_SEQ_FOR_EACH(YADE_HP_PY_COMPLEX_LOGGER, name, YADE_MINIEIGEN_HP) // instatinate templates for name<1>, name<2>, etc …
YADE_HP_PYTHON_COMPLEX_TEMPLATE_CREATE_LOGGER(ComplexVisitor);
