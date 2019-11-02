/**
 * \file unscented_kalman_filter.cc
 * \brief implementation of unscented_kalman_filter.h
*/

#include "unscented_kalman_filter.h"
#include "filter_model.h"

#include <iostream>
#include <cassert>
#include <stdexcept>

#include <limits>
#include <algorithm>

#include "QuadProg++_Eigen.h"
using namespace QuadProgPP;

UnscentedKalmanFilter::UnscentedKalmanFilter(SignalModel *filter_model)
    : m_FilterModel(filter_model), m_SigmaPointSpread(0.01)
{
  const unsigned int dim = m_FilterModel->state_dim();

  m_Scale = sqrt(dim + m_SigmaPointSpread);
  m_Weights(0) = m_SigmaPointSpread / (dim + m_SigmaPointSpread);
  // Create diagonal matrix.
  m_WeightsRepeated.setConstant(ukfZero);
  m_WeightsRepeated(0, 0) = m_Weights[0];
  for (unsigned int i = 1; i < (2 * dim) + 1; ++i)
  {
    m_Weights(i) = ukfHalf / (dim + m_SigmaPointSpread);
    m_WeightsRepeated(i, i) = m_Weights[i];
  }

  assert(static_cast<unsigned int>((m_FilterModel->Q()).rows()) == dim &&
         static_cast<unsigned int>((m_FilterModel->Q()).cols()) == dim);
}

void UnscentedKalmanFilter::SigmaPoints(const ukfStateVector &x,
                                        const ukfStateSquareMatrix &p,
                                        ukfStateCovMatrix &x_spread)
{
  // Horizontally stack x to the X_tmp matrix.
  ukfStateSquareMatrix X_tmp = x.rowwise().replicate(25); // CB: X changed to X_tmp to avoid notation confusion with member var X

  Eigen::LLT<ukfMatrixType> lltOfA(p);     // compute the Cholesky decomposition of A
  ukfStateSquareMatrix NewM = (lltOfA.matrixL()); // retrieve factor L  in the decomposition
  NewM *= m_Scale;
  // Create dim x (2 * dim + 1) matrix (x, x + m, x - m).
  x_spread.col(0) = x;
  x_spread.block<25, 25>(0, 1) = X_tmp + NewM;
  x_spread.block<25, 25>(0, 26) = X_tmp - NewM;
}

// vector version
void UnscentedKalmanFilter::Constrain(ukfStateVector &x, const ukfMatrixType &W)
{
  if (violatesContraints(x))
  {
    const ukfStateSquareMatrix WTranspose = W.transpose();
    ukfStateSquareMatrix W_tmp = (W + WTranspose) * ukfHalf;
    ukfStateVector g0 = -ukfOne * (W_tmp.transpose()) * x;
    const QPInequalityConstVec d = m_FilterModel->d(); // the inequality constraints
    const QPInequalityConst D = m_FilterModel->D(); // -- " --

    const ukfPrecisionType e = m_FilterModel->e(); // the equality constraints
    const ukfStateVector E = m_FilterModel->E(); // -- " --

    const ukfPrecisionType error = solve_quadprog(W_tmp, g0, E, e, D, d, x);
    //std::cout << "after " << x(21) << " " << x(22) << " " << x(23) << " sum " << x(21) + x(22) + x(23) << std::endl;
    if (error > 0.01) // error usually much smaller than that, if solve_quadprog fails it returns inf
    {
      throw std::logic_error("solve_quadprog error exceeds threshold 0.01!");
    }
  }
}

// matrix version
void UnscentedKalmanFilter::Constrain(ukfStateCovMatrix &localX, const ukfStateSquareMatrix &localW)
{
  for (unsigned int i = 0; i < 51; ++i)
  {
    ukfStateVector x = localX.col(i);
    Constrain(x, localW);
    localX.col(i) = x;
  }
}

bool UnscentedKalmanFilter::violatesContraints(ukfStateVector &x)
{
  const ukfVectorType d_test = (-ukfOne) * (m_FilterModel->D().transpose()) * x; // -D'*x
  for (unsigned int i = 0; i < d_test.size(); ++i)                               // if any(-D'*x > d) constraint is
                                                                                 // broken
  {
    if (d_test[i] > (m_FilterModel->d())[i])
    {
      return true;
    }
  }

  // Optimized for this UKF tractography only. Not a general implementation
  const ukfPrecisionType e_test = ((-ukfOne) * (m_FilterModel->E().transpose()) * x)[0]; // -E'*x
                                                                                         // if any(-E'*x != e) constraint is
                                                                                         // broken

  if (!(std::fabs(e_test - (m_FilterModel->e())) < std::numeric_limits<double>::epsilon()))
  {
    return true;
  }
  return false;
}

void UnscentedKalmanFilter::Filter(const ukfStateVector &x,
                                   const ukfStateSquareMatrix &p,
                                   const ukfVectorType &z,
                                   ukfStateVector &x_new,
                                   ukfStateSquareMatrix &p_new,
                                   ukfPrecisionType &dNormMSE)
{
  // Force a const version of the m_FilterModel to be used to ensure that it is not modified.
  SignalModel const *const localConstFilterModel = m_FilterModel;

  assert(static_cast<int>(x.size()) == localConstFilterModel->state_dim());
  assert(static_cast<int>(z.size()) == localConstFilterModel->signal_dim());
  assert(static_cast<int>(x_new.size()) == localConstFilterModel->state_dim());
  assert(static_cast<int>(p.rows()) == localConstFilterModel->state_dim() &&
         static_cast<int>(p.cols()) == localConstFilterModel->state_dim());
  assert(static_cast<int>(p_new.rows()) == localConstFilterModel->state_dim() &&
         static_cast<int>(p_new.cols()) == localConstFilterModel->state_dim());
  assert(localConstFilterModel);
  assert(static_cast<int>(m_Weights.size()) == 2 * localConstFilterModel->state_dim() + 1);

  const int signal_dim = localConstFilterModel->signal_dim();

  /** The state spread out according to the unscented transform */

  ukfStateCovMatrix X;
  X.setConstant(ukfZero);
  // Create sigma points.
  SigmaPoints(x, p, X); // doesnt change p, its const

  // X contains values > 1 for weights after sigma points
  if (localConstFilterModel->isConstrained())
  {
    // ukfMatrixType p_tmp = p; // will be changed in QuadProg
    Constrain(X, p);
  }

  // copy step needed because z is a const variable and can't be referenced
  ukfVectorType z_Eigen(z.size());
  for (unsigned int i = 0; i < z.size(); ++i)
  {
    z_Eigen[i] = z[i];
  }

  {
    if (localConstFilterModel->isRidgelets())
      localConstFilterModel->F(X, z_Eigen, p); // slightly negative fw is fixed here
    else
      localConstFilterModel->F(X); // slightly negative fw is fixed here
  }

  // std::cout<<"recorded signal: "<<z_Eigen<<std::endl;

  /** Used for the estimation of the new state */
  // std::cout << "\n X:"<<X <<"\n Weigths" << this->m_Weights;
  const ukfStateVector X_hat = X * this->m_Weights;
  ukfStateCovMatrix dim_dimext = X_hat.rowwise().replicate(51);

  const ukfStateSquareMatrix &Q = localConstFilterModel->Q();

  const ukfStateCovMatrix &X_ = X - dim_dimext;
  // Use const reference to avoid copying
  p_new = X_ * m_WeightsRepeated * X_.transpose() + Q;
  // std::cout<<"\n intitial P:"<<p_new;
  /** The signal */
  const ukfStateSquareMatrix Yk = p_new.inverse();
  const ukfStateVector yhat = Yk * X_hat;
  // std::cout << "\n P: " << p_new << "\n yhat: " << yhat <<"\n Q: "<<Q;

  ukfMatrixType Z(signal_dim, 51);
  Z.setConstant(ukfZero);

  localConstFilterModel->H(X, Z);

  /** Used for the estimation of the signal */
  ukfMatrixType signaldim_dimext(signal_dim, 51);
  signaldim_dimext.setConstant(ukfZero);
  const ukfVectorType Z_hat = Z * this->m_Weights;
  signaldim_dimext = Z_hat.rowwise().replicate(51);

  Z -= signaldim_dimext;
  const ukfMatrixType Z_ = Z;
  //const ukfMatrixType WeightsRepeated_Y_Transpose =m_WeightsRepeated*Y_.transpose();

  const ukfMatrixType temp = localConstFilterModel->R();
  ukfPrecisionType R = 1 / temp(0, 0);

  /** Covariance of the signal */
  //const ukfMatrixType Pyy = Y_ * WeightsRepeated_Y_Transpose + R;

  // Predict cross-correlation between state and observation.
  /** Covariance matrix state/signal */
  const ukfMatrixType Pxz = X_ * m_WeightsRepeated * Z_.transpose();

  // Kalman gain KalmanGainMatrix, estimate state/observation, compute covariance.
  // Solve KalmanGainMatrix = Pyy \ Pxy'
  // Solve Ax = b. Result stored in x. Matlab: x = A \ b.
  //x = A.ldlt().solve(b));
  const ukfMatrixType Ht = Yk * Pxz;

  const ukfVectorType zDiff = z_Eigen - Z_hat;
  dNormMSE = zDiff.squaredNorm() / z_Eigen.squaredNorm();

  const ukfMatrixType HtR = Ht * R;
  const ukfStateSquareMatrix I = HtR * Ht.transpose();
  const ukfStateVector i = HtR * (zDiff + (Pxz.transpose() * yhat));

  const ukfStateSquareMatrix YkI = Yk + I;
  p_new = YkI.inverse();
  // std::cout <<"\n p_new:"<<p_new;
  ukfStateVector x_new_Eigen = p_new * (i + yhat);

  if (localConstFilterModel->isConstrained())
  {
    Constrain(x_new_Eigen, YkI);
  }

  x_new = x_new_Eigen;
}
