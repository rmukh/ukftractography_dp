#ifndef RIDGELETS_INC_H_
#define RIDGELETS_INC_H_

#include "ukf_types.h"
#include "ukf_exports.h"

// Some ridgelets parameters hardcoded for now (need to add them in CLI?)
ukfPrecisionType sph_rho = 3.125;
unsigned int sph_J = 2;
ukfPrecisionType fista_lambda = 0.01;
unsigned int lvl = 4;
ukfPrecisionType max_odf_thresh = 0.7;

#endif