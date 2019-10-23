/**
 * \file fiber.h
 * \brief Description of a fiber
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#ifndef UKFFIBER_H_
#define UKFFIBER_H_

#include <vector>
#include <cassert>
#include "unscented_kalman_filter.h"
#include "linalg.h"

/**
 * \struct UKFFiber
 * \brief Points of a fiber, and scalars corresponding to the points
 *
 * The points that make a fiber are defined as a vector of 3D points. In addition
 * there is a vector of scalar values for each scalar value that can be recorded
 * of the same length
*/
struct UKFFiber
{
  /** vector of 3D points defining the fiber path */
  stdVec_t position;
  /** FA of tensor 1 */
  std::vector<ukfPrecisionType> fa;
  /** FA of tensor 2 */
  std::vector<ukfPrecisionType> fa2;
  /** FA of tensor 3 */
  std::vector<ukfPrecisionType> fa3;
  /** Array 2 norm of the covariance matrix */
  std::vector<ukfPrecisionType> norm;
  /** ukfStateVector of the current model at the current position*/
  std::vector<ukfStateVector> state;
  /** dim(state) x dim(state) matrix */
  std::vector<ukfMatrixType> covariance;
  /** Percentage of free water i.e. 1-w */
  std::vector<ukfPrecisionType> free_water;
  /** Normalized mean squared error of the signal reconstruction to the signal */
  std::vector<ukfPrecisionType> normMSE;
  /** Trace of tensor 1 */
  std::vector<ukfPrecisionType> trace;
  /** Trace of tensor 2 */
  std::vector<ukfPrecisionType> trace2;
  /** Weights (diffussion propagator, bi-exponential) */
  std::vector<ukfPrecisionType> w1;
  std::vector<ukfPrecisionType> w2;
  std::vector<ukfPrecisionType> w3;
  /** Angles (diffussion propagator, bi-exponential) */
  std::vector<ukfPrecisionType> w1w2angle;
  std::vector<ukfPrecisionType> w1w3angle;
  /** Uncertanity characteristics */
  std::vector<ukfPrecisionType> Fm1;
  std::vector<ukfPrecisionType> lmd1;
  std::vector<ukfPrecisionType> Fm2;
  std::vector<ukfPrecisionType> lmd2;
  std::vector<ukfPrecisionType> Fm3;
  std::vector<ukfPrecisionType> lmd3;
  std::vector<ukfPrecisionType> varW1;
  std::vector<ukfPrecisionType> varW2;
  std::vector<ukfPrecisionType> varW3;
  std::vector<ukfPrecisionType> varWiso;
};

/**
 * \brief Joins two fibers originating from the same seed point
 *
 * A pair of two primary fibers are started from each seed point in two opposite directions. This functions joins them up pairly to
 * form complete primary fibers, and eliminates fibers that are too short. Besides, each branch is back traced to form a whole fiber
*/
void PostProcessFibers(const std::vector<UKFFiber> &raw_primary, const std::vector<unsigned char> &discarded_fibers, std::vector<UKFFiber> &fibers);

/** The minimum number of points on a fiber. UKFFiber with fewer points are rejected */
const int MINIMUM_NUM_POINTS_ON_FIBER = 10;

/** Used to get rid of branches, that start near the end of primary fibers. See fiber.cc:70. */
const int FIBER_TAIL_THRESHOLD = 5;

#endif
