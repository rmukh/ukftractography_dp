/**
 * \file seed.h
 * \brief Contains data structure for storing seed point information.
*/

#ifndef SEED_H_
#define SEED_H_

#include "unscented_kalman_filter.h"
#include "linalg.h"

/**
 * \class SeedPointInfo
 * \brief Describes a seed point
 *
 * Stores all information for a seed point to start tractography. The
 * start_dir is only used for the simple model and the angles only for the
 * complex/full representation.
*/
class SeedPointInfo
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /** The state of the state-space represenation of the model. */
  ukfStateVector state;
  /** The covariance matrix of the state */
  ukfStateSquareMatrix covariance;
  /** The location of the seed */
  vec3_t point;
  /** The starting direction for the simple model */
  vec3_t start_dir;
  /** RTOP of the first compartment */
  ukfPrecisionType rtop1;
  /** RTOP of the second compartment */
  ukfPrecisionType rtop2;
  /** RTOP of the third compartment */
  ukfPrecisionType rtop3;
  /** RTOP from state */
  ukfPrecisionType rtop_model;
  /** RTOP from dMRI image */
  ukfPrecisionType rtop_signal;
};

/** Writes debug information about seeds to stdout. */
void PrintSeedInfo(const std::vector<SeedPointInfo> &seed_infos);

#endif
