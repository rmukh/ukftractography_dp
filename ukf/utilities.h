/**
 * \file utilities.h
 * \brief Calculation of frequently used (fa,ga,curve_radius) defined in this file
 * Calculations for NODDI model are implemented
*/
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>

#include "linalg.h"

#define NaN std::numeric_limits<ukfPrecisionType>::quiet_NaN()

/** Calculate fractional anisotropy from eigenvalues */
ukfPrecisionType l2fa(ukfPrecisionType l1, ukfPrecisionType l2, ukfPrecisionType l3);

/** Calculate Generalized anisotropy from signal */
ukfPrecisionType s2ga(const ukfMatrixType& signal);

/** Calculate mean signal function from signal */
ukfPrecisionType s2adc(const ukfMatrixType& signal);

/** Calculate curve radius from fiber */
ukfPrecisionType curve_radius(const stdVec_t& fiber);

#endif  // UTILITIES_H_
