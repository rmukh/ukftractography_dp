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
double Legendre(int n, double t)
{
  int k;
  double Pk_1, Pk_2, Pk; // P_{k-1}(x), P_{k-2}(x), P_k(x)

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

// n is half of maximum order of legendre polynomial by
void legendreGaussianIntegral(ukfVectorType ParComp, ukfPrecisionType n, ukfMatrixType &legGauInt)
{
  assert(n <= 6); // maximum value of n is 6

  ukfPrecisionType sqrtx, dx, emx;
  ukfMatrixType I(ParComp.size(), int(n + 1));
  ukfVectorType x2(ParComp.size());
  ukfVectorType x3(ParComp.size());
  ukfVectorType x4(ParComp.size());
  ukfVectorType x5(ParComp.size());
  ukfVectorType x6(ParComp.size());
  for (int i = 0; i < ParComp.size(); i++)
  {
    if (ParComp[i] > 0.05)
    {
      sqrtx = sqrt(ParComp[i]);
      I(i, 0) = sqrt(M_PI) * erf(sqrtx) / sqrtx;
      dx = 1 / ParComp[i];
      emx = -1 * exp(-1 * ParComp[i]);
      for (int j = 1; j < n + 1; j++)
      {
        I(i, j) = emx + (j - 0.5) * I(i, j - 1);
        I(i, j) = I(i, j) * dx;
      }
    }
    else
    {
      I(i, 0) = 0;
      x2[i] = ParComp[i] * ParComp[i];
      x3[i] = x2[i] * ParComp[i];
      x4[i] = x3[i] * ParComp[i];
      x5[i] = x4[i] * ParComp[i];
      x6[i] = x5[i] * ParComp[i];
    }
  }
  for (int i = 0; i < ParComp.size(); i++)
  {
    if (ParComp[i] > 0.05)
    {
      for (int j = 0; j < n + 1; j++)
      {
        if (j == 0)
          legGauInt(i, 0) = I(i, 0);
        else if (j == 1)
          legGauInt(i, 1) = -0.5 * I(i, 0) + 1.5 * I(i, 1);
        else if (j == 2)
          legGauInt(i, 2) = 0.375 * I(i, 0) - 3.75 * I(i, 1) + 4.375 * I(i, 2);
        else if (j == 3)
          legGauInt(i, 3) = -0.3125 * I(i, 0) + 6.5625 * I(i, 1) - 19.6875 * I(i, 2) + 14.4375 * I(i, 3);
        else if (j == 4)
          legGauInt(i, 4) = 0.2734375 * I(i, 0) - 9.84375 * I(i, 1) + 54.140625 * I(i, 2) - 93.84375 * I(i, 3) + 50.2734375 * I(i, 4);
        else if (j == 5)
          legGauInt(i, 5) = -0.24609 * I(i, 0) + 13.53515 * I(i, 1) - 117.30468 * I(i, 2) + 351.91406 * I(i, 3) - 427.32421 * I(i, 4) + 180.425781 * I(i, 5);
        else if (j == 6)
          legGauInt(i, 6) = 0.22558 * I(i, 0) - 17.595703 * I(i, 1) + 219.946289 * I(i, 2) - 997.089843 * I(i, 3) + 2029.79003 * I(i, 4) - 1894.470703 * I(i, 5) + 660.194335 * I(i, 6);
      }
    }
    else
    {
      for (int j = 0; j < n + 1; j++)
      {
        if (j == 0)
          legGauInt(i, 0) = 2 - 2 * ParComp[i] / 3 + x2[i] / 5 - x3[i] / 21 + x4[i] / 108;
        else if (j == 1)
          legGauInt(i, 1) = -4 * ParComp[i] / 15 + 4 * x2[i] / 35 - 2 * x3[i] / 63 + 2 * x4[i] / 297;
        else if (j == 2)
          legGauInt(i, 2) = 8 * x2[i] / 315 - 8 * x3[i] / 693 + 4 * x4[i] / 1287;
        else if (j == 3)
          legGauInt(i, 3) = -16 * x3[i] / 9009 + 16 * x4[i] / 19305;
        else if (j == 4)
          legGauInt(i, 4) = 32 * x4[i] / 328185;
        else if (j == 5)
          legGauInt(i, 5) = -64 * x5[i] / 14549535;
        else if (j == 6)
          legGauInt(i, 6) = 128 * x6[i] / 760543875;
      }
    }
  }
}

void ParallelComponent(ukfPrecisionType dPar,
                       ukfVectorType &gradientStrength,
                       ukfVectorType &pulseSeparation,
                       ukfVectorType &ParComp)
{
  ukfPrecisionType GAMMA, modQ_sq, difftime;
  GAMMA = 267598700;

  for (int i = 0; i < gradientStrength.size(); i++)
  {
    modQ_sq = GAMMA * GAMMA * pulseSeparation[i] * pulseSeparation[i] * gradientStrength[i] * gradientStrength[i];
    difftime = pulseSeparation[i] - pulseSeparation[i] / 3;
    ParComp[i] = modQ_sq * difftime * dPar;
  }
}

