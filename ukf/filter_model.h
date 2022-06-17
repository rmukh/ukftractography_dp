/**
 * \struct RIDG
 * \brief Ridgelets model
 *
 * Model describing tractography with ridgelets with free water
 */

#ifndef FILTER_MODEL_H_
#define FILTER_MODEL_H_

#include <cassert>
#include <vector>
#include <iostream>

#include "linalg.h"
#include "ISignalData.h"
#include "unscented_kalman_filter.h"

#include "math_utilities.h"

#include "SOLVERS.h"
#include "SPH_RIDG.h"
#include "UtilMath.h"

inline mat33_t SetIdentityScaled(ukfPrecisionType diff_fw)
{
  mat33_t tmp;
  tmp.setIdentity();
  tmp *= diff_fw;
  return tmp;
}

const ukfPrecisionType Pi(std::atan(static_cast<ukfPrecisionType>(1.0)) * 4);

inline ukfPrecisionType BhattacharyyaCoeff(ukfPrecisionType x_sr, ukfPrecisionType x_pred, ukfPrecisionType cov)
{
  return std::exp(-0.25 * (std::pow((x_sr - x_pred), 2) / cov));
}

inline ukfPrecisionType BhattacharyyaCoeff(vec3_t &x_sr, vec3_t &x_pred, const ukfMatrixType &cov)
{
  vec3_t diff = x_sr - x_pred;
  return std::exp(-0.125 * diff.transpose() * cov.inverse() * diff);
}

inline ukfPrecisionType BhattacharyyaCoeff(vec3_t &x_sr, vec3_t &x_pred, const ukfMatrixType &cov, const ukfMatrixType &cov2)
{
  vec3_t diff = x_sr - x_pred;
  ukfMatrixType cov_total = (cov + cov2) / 2.0;
  return std::exp(-((0.125 * diff.transpose() * cov_total.inverse() * diff) + 0.5 * std::log(cov_total.determinant() / std::sqrt(cov.determinant() * cov2.determinant()))));
}

inline ukfPrecisionType AngularSimilarity(vec3_t &x_sr, vec3_t &x_pred)
{
  ukfPrecisionType dot = x_sr.dot(x_pred);
  ukfPrecisionType den_a = x_sr.norm();
  ukfPrecisionType den_b = x_pred.norm();

  if (cmpf(den_a, ukfZero) || cmpf(den_b, ukfZero))
  {
    std::cout << "cosine similarity is not defined whenever one or both input vectors are zero-vectors." << std::endl;
    throw;
  }

  return ukfOne - (std::acos(std::min(std::max(dot / (den_a * den_b), -ukfOne), ukfOne)) / Pi);
}

/**
 * \struct SignalModel
 * \brief Implementation of a signal model
 * A class that defines the transition function and the observation
 * model to be used in the Kalman filter.
 */
class SignalModel
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** Constructor */
  SignalModel(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qt, ukfPrecisionType qw, ukfPrecisionType qwiso,
              ukfPrecisionType rs, bool constrained, const ukfPrecisionType diff_fw,
              ukfMatrixType &Aridg, ukfMatrixType &Qridg, ukfMatrixType &fcsridg, ukfMatrixType &nuridg,
              vector<vector<unsigned>> &connridg, signalMaskType &sm, ukfPrecisionType fl, ukfPrecisionType mot)
      : _state_dim(25), _rs(rs), _signal_dim(0), _signal_data(NULL), _constrained(constrained), _ridgelets_used(false),
        _lambda_min_fast_diffusion(1.0), _lambda_min_slow_diffusion(0.1), _lambda_max_diffusion(3000.0),
        _w_fast_diffusion(0.7), m_D_iso(SetIdentityScaled(diff_fw)), A(Aridg), QRidg(Qridg), fcs(fcsridg), nu(nuridg),
        conn(connridg), signal_mask(sm), fista_lambda(fl), max_odf_thresh(mot)
  {
    // size(X, 'column') == 25
    // X = [x10, x11, x12, l11, l12, l13, l14, x20, x21, x22, l21, l22, l23, l24, x30, x31, x32, l31, l32, l33, l34, w1, w2, w3, wiso]'
    //     [0  , 1  , 2  , 3  , 4  , 5  , 6  , 7  , 8  , 9  , 10 , 11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 , 22, 23, 24]
    // There are three tensors at max in this model, and we have a bi-exponential model

    // Q is the injected process noise bias. It is a diagonal matrix. It is used in the state transfer function.

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(7, 7) = _Q(8, 8) = _Q(9, 9) = _Q(14, 14) = _Q(15, 15) = _Q(16, 16) = qs; // noise for the main direction
    _Q(3, 3) = _Q(4, 4) = _Q(10, 10) = _Q(11, 11) = _Q(17, 17) = _Q(18, 18) = ql;                                // noise for the lambdas (fast diffusion)
    _Q(5, 5) = _Q(6, 6) = _Q(12, 12) = _Q(13, 13) = _Q(19, 19) = _Q(20, 20) = qt;                                // noise for the lambdas (slow diffusion)
    _Q(21, 21) = _Q(22, 22) = _Q(23, 23) = qw;                                                                   // noise for the omega's
    _Q(24, 24) = qwiso;                                                                                          // noise for the free water weight

    // D is the constraint matrix.
    // The 1st dim of D is used to select which component of x we want to constraint, the 2nd is for the inequality (*-1 -> reverse the inequality)
    // d is the contraint value
    // D'*x >= -d

    // N_constr constraints for the 25 dimensions of the state
    _D.setConstant(ukfZero);

    /* Setting the constraints according to D'*x >= -d */

    // Tensor 1 (minimal values)
    _D(3, 0) = 1;
    _d(0) = _lambda_min_fast_diffusion; // l11 >= min
    _D(4, 1) = 1;
    _d(1) = _lambda_min_fast_diffusion; // l12 >= min
    _D(5, 2) = 1;
    _d(2) = _lambda_min_slow_diffusion; // l13 >= min
    _D(6, 3) = 1;
    _d(3) = _lambda_min_slow_diffusion; // l14 >= min

    // Tensor 2 (minimal values)
    _D(10, 4) = 1;
    _d(4) = _lambda_min_fast_diffusion; // l21 >= min
    _D(11, 5) = 1;
    _d(5) = _lambda_min_fast_diffusion; // l22 >= min
    _D(12, 6) = 1;
    _d(6) = _lambda_min_slow_diffusion; // l23 >= min
    _D(13, 7) = 1;
    _d(7) = _lambda_min_slow_diffusion; // l24 >= min

    // Tensor 3 (minimal values)
    _D(17, 8) = 1;
    _d(8) = _lambda_min_fast_diffusion; // l31 >= min
    _D(18, 9) = 1;
    _d(9) = _lambda_min_fast_diffusion; // l32 >= min
    _D(19, 10) = 1;
    _d(10) = _lambda_min_slow_diffusion; // l33 >= min
    _D(20, 11) = 1;
    _d(11) = _lambda_min_slow_diffusion; // l34 >= min

    // Tensor 1 (maximal values)
    _D(3, 12) = -1; // l11 <= max
    _d(12) = _lambda_max_diffusion;
    _D(4, 13) = -1; // l12 <= max
    _d(13) = _lambda_max_diffusion;
    _D(5, 14) = -1; // l13 <= max
    _d(14) = _lambda_max_diffusion;
    _D(6, 15) = -1; // l14 <= max
    _d(15) = _lambda_max_diffusion;

    // Tensor 2 (maximal values)
    _D(10, 16) = -1; // l21 <= max
    _d(16) = _lambda_max_diffusion;
    _D(11, 17) = -1; // l22 <= max
    _d(17) = _lambda_max_diffusion;
    _D(12, 18) = -1; // l23 <= max
    _d(18) = _lambda_max_diffusion;
    _D(13, 19) = -1; // l24 <= max
    _d(19) = _lambda_max_diffusion;

    // Tensor 3 (maximal values)
    _D(17, 20) = -1;
    _d(20) = _lambda_max_diffusion; // l31 <= max
    _D(18, 21) = -1;
    _d(21) = _lambda_max_diffusion; // l32 <= max
    _D(19, 22) = -1;
    _d(22) = _lambda_max_diffusion; // l33 <= max
    _D(20, 23) = -1;
    _d(23) = _lambda_max_diffusion; // l34 <= max

    // Weights
    _D(21, 24) = 1;
    _d(24) = 0; // w1 >= 0
    _D(21, 25) = -1;
    _d(25) = 1; // w1 <= 1

    _D(22, 26) = 1;
    _d(26) = 0; // w2 >= 0
    _D(22, 27) = -1;
    _d(27) = 1; // w2 <= 1

    _D(23, 28) = 1;
    _d(28) = 0; // w3 >= 0
    _D(23, 29) = -1;
    _d(29) = 1; // w3 <= 1

    // Free water
    _D(24, 30) = 1;
    _d(30) = 0; // wiso >= 0
    _D(24, 31) = -1;
    _d(31) = 1; // wiso <= 1

    // Equality constraints (w1 + w2 + w3 = 1)
    _E.setConstant(ukfZero);

    _E(21, 0) = _E(22, 0) = _E(23, 0) = -1.0;
    _e = 1.0;
  }

  /** Destructor */
  ~SignalModel()
  {
  }

  /** state transition function */
  void F(ukfStateCovMatrix & /* X */, ukfVectorType /* s */, const ukfMatrixType & /* &covMatrix */) const;
  void F(ukfStateCovMatrix & /* X */) const;

  /** observation, i.e. signal reconstruction */
  void H(const ukfStateCovMatrix & /* X */, ukfMatrixType & /* Y */) const;
  void H(const ukfStateVector &, ukfMatrixType &) const;

  /* 
  Functions that convert the state into a tensor representation that can be
  used for tractography (meaning the main direction and all the eigenvalues/lambdas).
  */

  /** Extracts principal diffusion direction and eigen values from the state */
  void State2Tensor3T(const ukfStateVector &x, const vec3_t &old_m, vec3_t &m1, vec3_t &l1, vec3_t &m2, vec3_t &l2, vec3_t &m3, vec3_t &l3);
  void State2Tensor3T(const ukfStateVector &x, const vec3_t &old_m, vec3_t &m1);

  ukfPrecisionType cosine_similarity(vec3_t &First, vec3_t &Second) const;
  /** The minimum/maximum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min_fast_diffusion;
  const ukfPrecisionType _lambda_min_slow_diffusion;
  const ukfPrecisionType _lambda_max_diffusion;

  /** The weight fractions of the fast diffusion components */
  const ukfPrecisionType _w_fast_diffusion;

  /** apparent diffusion coefficient of free water */
  mat33_t m_D_iso;

  /** Ridgelets matricies/vectors **/
  ukfMatrixType &A;
  ukfMatrixType &QRidg;
  ukfMatrixType &fcs;
  ukfMatrixType &nu;
  vector<vector<unsigned>> &conn;

  signalMaskType &signal_mask;

  const ukfPrecisionType fista_lambda;
  const ukfPrecisionType max_odf_thresh;

  /** Returns the dimension of the state */
  int state_dim() const
  {
    return _state_dim;
  }

  /** Set the dimension of the signal, and resize the corresponding matrices correspondigly */
  void set_signal_dim(const int dim)
  {
    _signal_dim = dim;

    _R.resize(_signal_dim, _signal_dim);
    _R.setConstant(ukfZero); // necessary because otherwise there is memory leftovers in the matrix
    for (int i = 0; i < _signal_dim; ++i)
    {
      _R(i, i) = _rs;
    }
  }

  /** Returns the dimension of the signal */
  int signal_dim() const
  {
    return _signal_dim;
  }

  /** The noise in the state transfer function used by the UKF. */
  const ukfStateSquareMatrix &Q() const
  {
    return _Q;
  }

  /** The noise in the signal reconstruction used by the UKF. */
  const ukfMatrixType &R() const
  {
    return _R;
  }

  // The next two functions are only used for the constrained case

  /** The inequality constraint matrix for the constrained UKF */
  const QPInequalityConst &D() const
  {
    return _D;
  }

  /** The inequality constraint right hand side for the constrained UKF */
  const QPInequalityConstVec &d() const
  {
    return _d;
  }

  /** The equality constraint matrix for the constrained UKF */
  const ukfStateVector &E() const
  {
    return _E;
  }

  /** The equality constraint right hand side for the constrained UKF */
  const ukfPrecisionType &e() const
  {
    return _e;
  }

  /** Are we running the contrained version of the filter? */
  bool isConstrained() const
  {
    return _constrained;
  }

  /** Set the pointer to the diffusion signal data */
  void set_signal_data(ISignalData *signal_data)
  {
    _signal_data = signal_data;
  }

  /** Are we using Spherical Ridgelets **/
  bool isRidgelets() const
  {
    return _ridgelets_used;
  }

protected:
  /** Checks if d is smaller than a small negative threshold. If yes an error is returned. Otherwise d is rounded to ukfZero
   */
  ukfPrecisionType CheckZero(const ukfPrecisionType &d, const std::string &func_name) const;

  /** The dimension of the state */
  const int _state_dim;

  /** The constant signal noise for each signal component */
  ukfPrecisionType _rs;

  /** The dimension of the signal */
  int _signal_dim;

  /** Pointer to the diffusion signal data */
  ISignalData *_signal_data;

  /** Process noise */
  ukfStateSquareMatrix _Q;

  /** Signal reconstruction noise */
  ukfMatrixType _R;

  /** Inequality constraint matrix, only used for constrained UKF */
  QPInequalityConst _D;

  /** Inequality right hand side, only used for constrained UKF */
  QPInequalityConstVec _d;

  /** Equality constraint matrix, only used for constrained UKF */
  ukfStateVector _E;

  /** Equality right hand side, only used for constrained UKF */
  ukfPrecisionType _e;

  /** Are we using the constrained filter */
  bool _constrained;

  /** Are we using spherical ridgelets**/
  bool _ridgelets_used;
};

#endif // FILTER_MODEL_H_
