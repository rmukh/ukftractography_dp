/**
 * \file utilities.cc
 * \brief implementation of utilities.h
*/

#include "utilities.h"
#include "math_utilities.h"

// std includes
#include <math.h>
#include <float.h>

ukfPrecisionType l2fa(ukfPrecisionType l1, ukfPrecisionType l2, ukfPrecisionType l3)
{
  if (std::fabs(l2 - l3) < std::numeric_limits<ukfPrecisionType>::epsilon())
  {
    return std::fabs(l1 - l2) / sqrt(l1 * l1 + 2.0 * l2 * l2);
  }
  else
  {
    return sqrt(ukfHalf * ((l1 - l2) * (l1 - l2) + (l2 - l3) * (l2 - l3) + (l1 - l3) * (l3 - l1)) / (l1 * l1 + l2 * l2 + l3 * l3));
  }
}

ukfPrecisionType s2ga(const ukfMatrixType &signal)
{

  int n = static_cast<int>(signal.rows());

  assert(signal.cols() == 1);

  ukfPrecisionType mu = ukfZero;
  for (int i = 0; i < n; ++i)
  {
    mu += signal(i, 0);
  }
  // average
  mu = mu / n;

  ukfPrecisionType mu_sq = ukfZero;
  ukfPrecisionType mu_sub = ukfZero;
  for (int i = 0; i < n; ++i)
  {
    mu_sq += signal(i, 0) * signal(i, 0);
    mu_sub += (signal(i, 0) - mu) * (signal(i, 0) - mu);
  }

  return sqrt(mu_sub * n) / sqrt((n - 1) * mu_sq);
}

// mean signal function
ukfPrecisionType s2adc(const ukfMatrixType &signal)
{
  unsigned n = static_cast<int>(signal.rows());

  assert(signal.cols() == 1);

  ukfPrecisionType mu = ukfZero;
  for (unsigned i = 0; i < n; ++i)
  {
    mu += signal(i, 0);
  }
  // average
  mu = mu / n;
  return mu;
}

ukfPrecisionType curve_radius(const stdVec_t &fiber)
{
  size_t length = fiber.size();

  if (length < 3)
  {
    return ukfOne;
  }

  vec3_t v1 = fiber[length - 2] - fiber[length - 3];
  vec3_t v2 = fiber[length - 1] - fiber[length - 2];

  // Normalize
  v1.normalize();
  v2.normalize();
  ukfPrecisionType n1 = v1.norm();
  ukfPrecisionType n2 = v2.norm();

  ukfPrecisionType curv = ((v2 - v1) / (n2 + n1)).norm();
  if (std::isnan(curv))
  {
    return ukfZero;
  }

  return ukfOne / curv;
}

// legendre polynomial of degree n and order 0
ukfPrecisionType Legendre(int n, ukfPrecisionType t)
{
  int k;
  ukfPrecisionType Pk_1, Pk_2, Pk; // P_{k-1}(x), P_{k-2}(x), P_k(x)

  Pk_2 = 0.0;
  Pk_1 = 1.0;
  Pk = 1.0;

  for (k = 1; k <= n; k++)
  {
    Pk = (2.0 * k - 1.0) / k * t * Pk_1 - (k - 1.0) / k * Pk_2;
    Pk_2 = Pk_1;
    Pk_1 = Pk;
  }

  return Pk;
}

