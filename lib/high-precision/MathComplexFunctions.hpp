/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This file contains mathematical functions available in standard library and boost library.
//     https://en.cppreference.com/w/cpp/numeric/math
//     https://en.cppreference.com/w/cpp/numeric/special_functions
// They have to be provided here as inline redirections towards the correct implementation, depending on what precision type yade is being compiled with.
// This is the only way to make sure that ::std, ::boost::math, ::boost::multiprecision are all called correctly.

// TODO: Boost documentation recommends to link with tr1: -lboost_math_tr1 as it provides significant speedup. For example replace boost::math::acosh(x) ↔ boost::math::tr1::acosh(x)
//     https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/overview_tr1.html
//#include <boost/math/tr1.hpp>

// All functions have to be tested in:
// * py/high-precision/_math.cpp
// * py/tests/testMath.py
// * py/tests/testMathHelper.py

#ifndef YADE_MATH_COMPLEX_FUNCIONS_HPP
#define YADE_MATH_COMPLEX_FUNCIONS_HPP

#include <lib/high-precision/Real.hpp>

namespace forCtags {
struct MathComplexFunctions {
}; // for ctags
}

/*************************************************************************/
/*********************** YADE_REAL_MATH_NAMESPACE ************************/
/*************************************************************************/
// The YADE_REAL_MATH_NAMESPACE macro makes sure that proper namespace is used for Real and eventual RealHP<…> types.
// On older systems RealHP<…> can't work, so we need to check if it was explicitly disabled during cmake autodetection.
#if ((YADE_REAL_BIT <= 80) and defined(YADE_DISABLE_REAL_MULTI_HP))
#define YADE_REAL_MATH_NAMESPACE ::std
#else
#define YADE_REAL_MATH_NAMESPACE ::boost::multiprecision
#endif

namespace yade {

namespace complex_literals {

	// https://en.cppreference.com/w/cpp/numeric/complex/operator%22%22i , that is C++14 native operator ""i.
	// But operator "" literal without starting _ are reserved for C++ standard. Se we provide here ""_i.
	// And you can write:
	//   using namespace std::complex_literals;
	//   auto    a = 5i;        // and have a std::complex<double> 5i   , without '_' it's the natice C++14 literal.
	//   auto    b = 5.0i;      // and have a std::complex<double> 5.0i

	//   auto    a = 5_i;       // and have a      Complex         5_i  , with '_' it's our own custom literal.
	//   auto    b = 5.0_i;     // and have a      Complex         5.0_i
	//   Complex c = 5_i;       // and have a      Complex         5_i
	//   Complex d = 5.0_i;     // and have a      Complex         5.0_i
	//
	// NOTE: https://stackoverflow.com/questions/16596864/c11-operator-with-double-parameter
	//       Only two overloads:
	//         * long double
	//         * unsigned long long int
	//       Because these are literals. To be typed directly into the code: 5_i or 5.0_i. No ambiguity allowed.
	inline Complex operator""_i(const long double a) { return Complex { Real { 0 }, static_cast<Real>(a) }; }
	inline Complex operator""_i(const unsigned long long int a) { return Complex { Real { 0 }, static_cast<Real>(a) }; }

} // namespace complex_literals

using namespace complex_literals;

namespace math {

	/********************************************************************************************/
	/**********************          Complex conj, real, imag, abs          *********************/
	/********************************************************************************************/
	// Add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp and py/tests/testMath.py
	// Note: most of the functions above can be converted to accept complex by changing levelOfRealHP<> → levelOfHP<>, provided that a complex version exists.
	// The check involving int Level = levelOfComplexHP<Cc> is necessary to make sure that function is called only with yade Complex supported HP types.
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline Cc conj(const Cc& a)
	{
		using ::std::conj;
		using YADE_REAL_MATH_NAMESPACE::conj;
		return conj(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline RealHP<Level> real(const Cc& a)
	{
		using ::std::real;
		using YADE_REAL_MATH_NAMESPACE::real;
		return real(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline RealHP<Level> imag(const Cc& a)
	{
		using ::std::imag;
		using YADE_REAL_MATH_NAMESPACE::imag;
		return imag(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, RealHP<Level>>::type abs(const Cc& a)
	{ // # Warning: C++ std::norm is squared length. In python/mpmath norm is C++ std::abs == length.
		using ::std::abs;
		using YADE_REAL_MATH_NAMESPACE::abs;
		return abs(static_cast<const UnderlyingHP<Cc>&>(a));
	}

	/********************************************************************************************/
	/**********************          Complex arg, norm, proj, polar         *********************/
	/********************************************************************************************/
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline RealHP<Level> arg(const Cc& a)
	{
		using ::std::arg;
		using YADE_REAL_MATH_NAMESPACE::arg;
		return arg(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline RealHP<Level> norm(const Cc& a)
	{ // note: Returns the squared magnitude of the complex number z, see https://en.cppreference.com/w/cpp/numeric/complex/norm
		// # Warning: C++ norm is squared length. In python/mpmath norm is C++ abs == length.
		using ::std::norm;
		using YADE_REAL_MATH_NAMESPACE::norm;
		return norm(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline Cc proj(const Cc& a)
	{
		using ::std::proj;
		using YADE_REAL_MATH_NAMESPACE::proj;
		return proj(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Rr, int Level = levelOfRealHP<Rr>> inline ComplexHP<Level> polar(const Rr& rho, const Rr& theta)
	{ // https://en.cppreference.com/w/cpp/numeric/complex/polar
		using ::std::polar;
		using YADE_REAL_MATH_NAMESPACE::polar;
		// Use tmp to ensure proper ComplexHP return type. Usually it's just optimised away. Rarely it fixes a problem with boost::multiprecision using a different complex type.
		auto tmp = polar(static_cast<const UnderlyingHP<Rr>&>(rho), static_cast<const UnderlyingHP<Rr>&>(theta));
		return { tmp.real(), tmp.imag() };
	}

	/********************************************************************************************/
	/**********************         Complex trigonometric functions        **********************/
	/********************************************************************************************/
	// typename RC is a type which can be Real or Complex, Rr → only Real, Cc → only Complex.
	// The check involving int Level = levelOfComplexHP<Cc> is necessary to make sure that function is called only with yade supported HP types.
	// int Level is the N in RealHP<N>.
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type sin(const Cc& a)
	{
		using ::std::sin;
		using YADE_REAL_MATH_NAMESPACE::sin;
		return sin(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type sinh(const Cc& a)
	{
		using ::std::sinh;
		using YADE_REAL_MATH_NAMESPACE::sinh;
		return sinh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type cos(const Cc& a)
	{
		using ::std::cos;
		using YADE_REAL_MATH_NAMESPACE::cos;
		return cos(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type cosh(const Cc& a)
	{
		using ::std::cosh;
		using YADE_REAL_MATH_NAMESPACE::cosh;
		return cosh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type tan(const Cc& a)
	{ // Warning: in boost < 1.76 the float128 and higher complex tan(…), tanh(…) has huge ULP error.
		using ::std::tan;
		using YADE_REAL_MATH_NAMESPACE::tan;
		return tan(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type tanh(const Cc& a)
	{ // Warning: in boost < 1.76 the float128 and higher complex tan(…), tanh(…) has huge ULP error.
		using ::std::tanh;
		using YADE_REAL_MATH_NAMESPACE::tanh;
		return tanh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	// add more Real or Complex functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py

	/********************************************************************************************/
	/**********************     Complex inverse trigonometric functions    **********************/
	/********************************************************************************************/
	// The check involving int Level = levelOfComplexHP<Cc> is necessary to make sure that function is called only with yade Real supported HP types.
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type asin(const Cc& a)
	{
		using ::std::asin;
		using YADE_REAL_MATH_NAMESPACE::asin;
		return asin(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type asinh(const Cc& a)
	{
		using ::std::asinh;
		using YADE_REAL_MATH_NAMESPACE::asinh;
		return asinh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type acos(const Cc& a)
	{
		using ::std::acos;
		using YADE_REAL_MATH_NAMESPACE::acos;
		return acos(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type acosh(const Cc& a)
	{
		using ::std::acosh;
		using YADE_REAL_MATH_NAMESPACE::acosh;
		return acosh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type atan(const Cc& a)
	{
		using ::std::atan;
		using YADE_REAL_MATH_NAMESPACE::atan;
		return atan(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type atanh(const Cc& a)
	{
		using ::std::atanh;
		using YADE_REAL_MATH_NAMESPACE::atanh;
		return atanh(static_cast<const UnderlyingHP<Cc>&>(a));
	}

	/********************************************************************************************/
	/**********************   logarithm, exponential and power functions   **********************/
	/********************************************************************************************/
	// Add more functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	// They can be converted to accept complex by changing levelOfRealHP<> → levelOfHP<>, provided that a complex version exists.
	// But remember to add tests for complex versions in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type exp(const Cc& a)
	{
		using ::std::exp;
		using YADE_REAL_MATH_NAMESPACE::exp;
		return exp(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type log(const Cc& a)
	{
		using ::std::log;
		using YADE_REAL_MATH_NAMESPACE::log;
		return log(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type log10(const Cc& a)
	{
		using ::std::log10;
		using YADE_REAL_MATH_NAMESPACE::log10;
		return log10(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename A, typename B, int Level = levelOfComplexHP<A>, typename Cc = PromoteHP<A>>
	inline typename boost::enable_if_c<(std::is_convertible<B, Cc>::value and isComplexHP<A>), Cc>::type pow(const A& a, const B& b)
	{
		using ::std::pow;
		using YADE_REAL_MATH_NAMESPACE::pow;
		return pow(static_cast<const UnderlyingHP<Cc>&>(a), static_cast<const UnderlyingHP<Cc>&>(b));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type sqrt(const Cc& a)
	{
		using ::std::sqrt;
		using YADE_REAL_MATH_NAMESPACE::sqrt;
		return sqrt(static_cast<const UnderlyingHP<Cc>&>(a));
	}


} // namespace math
} // namespace yade

#undef YADE_REAL_MATH_NAMESPACE
#endif
