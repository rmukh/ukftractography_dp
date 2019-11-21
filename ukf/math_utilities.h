/**
* \file math_utilities.h
* \brief Implements some math functions needed in various parts of the code
*
* This file reimplements some functionality from math.h in order to overcome some cross-platform problems.
* There are no execution speed repercussions resulting from this substitution.
*
* \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef MATH_UTILITIES_H_
#define MATH_UTILITIES_H_

#include <limits>

#ifndef M_PI
// Source: http://stackoverflow.com/questions/13690483/better-more-portable-method-of-defining-pi-in-c-c
//         http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html

#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
#endif

#define UKF_PI M_PI
#define PI_COEFF 5.5683279968317078452848179821188357020136243902832439 // ~ std::pow(UKF_PI, 1.5);

#define DEG_TO_RAD (UKF_PI / 180.0)
#define RAD_TO_DEG (180.0 / UKF_PI)
inline ukfPrecisionType DegToRad(const ukfPrecisionType deg)
{
   return deg * DEG_TO_RAD;
}

inline ukfPrecisionType RadToDeg(const ukfPrecisionType rad)
{
   return rad * RAD_TO_DEG;
}

/* Relative error bounded by 3e-9 for normalized outputs
   Returns invalid outputs for nan inputs
   Continuous error
   Vectorizable only with AVX512dq extensions because of the
   double->int64 cast. On GCC, use option -mavx512dq. */
inline double expapprox_d(double val)
{
   double exp_cst1_d = 9218868437227405312.;
   double exp_cst2_d = 0.;

   union {
      int64_t i;
      double f;
   } xu, xu2;
   double val2, val3, val4, b;
   int64_t val4i;
   val2 = 6497320848556798.092 * val + 4607182418800017408.;
   val3 = val2 < exp_cst1_d ? val2 : exp_cst1_d;
   val4 = val3 > exp_cst2_d ? val3 : exp_cst2_d;
   val4i = (int64_t)val4;
   xu.i = val4i & 0x7FF0000000000000;
   xu2.i = (val4i & 0xFFFFFFFFFFFFF) | 0x3FF0000000000000;
   b = xu2.f;

   /* Generated in Sollya with:
     > f=remez(1-x*exp(-(x-1)*log(2)),
               [|(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
                 (x-1)*(x-2)*x*x*x, (x-1)*(x-2)*x*x*x*x|],
                 [1.000001,1.999999], exp(-(x-1)*log(2)));
     > plot(exp((x-1)*log(2))/(f+x)-1, [1,2]);
     > f+x;
  */
   return xu.f * (0.5002494548377929861 + b * (0.3453457447382168695 + b *
                                                                           (0.1226618159032020501 + b * (2.4869768911673215212e-2 + b *
                                                                                                                                        (6.7148791796145116483e-3 + b * (-5.8813825553864185693e-5 + b *
                                                                                                                                                                                                         2.17150255054231565039e-4))))));
}

#endif // MATH_UTILITIES_H_
