#ifndef RIDG_BiExp_FW__
#define RIDG_BiExp_FW__

#include "filter_model.h"

#include "SOLVERS.h"
#include "SPH_RIDG.h"
#include "UtilMath.h"
/**
 * \struct RIDG
 * \brief Ridgelets model
 *
 * Model describing tractography with ridgelets with free water
*/

class Ridg_BiExp_FW : public FilterModel
{
public:
    Ridg_BiExp_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qt, ukfPrecisionType qw, ukfPrecisionType qwiso,
                  ukfPrecisionType rs, const ukfVectorType &weights_on_tensors, bool constrained, const ukfPrecisionType diff_fw,
                  ukfMatrixType &Aridg, ukfMatrixType &Qridg, ukfMatrixType &fcsridg, ukfMatrixType &nuridg,
                  vector<vector<unsigned>> &connridg, signalMaskType &sm, ukfPrecisionType fl, ukfPrecisionType mot)
        : FilterModel(25, rs, weights_on_tensors, constrained, true),
          _lambda_min_fast_diffusion(1.0), _lambda_min_slow_diffusion(0.1), _lambda_max_diffusion(3000.0),
          _w_fast_diffusion(0.7), m_D_iso(SetIdentityScaled(diff_fw)), A(Aridg), Q(Qridg), fcs(fcsridg), nu(nuridg), conn(connridg),
          signal_mask(sm), fista_lambda(fl), max_odf_thresh(mot)
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

        const unsigned int N_constr = 32;

        // N_constr constraints for the 25 dimensions of the state
        _D.resize(25, N_constr);
        _D.setConstant(ukfZero);

        _d.resize(N_constr);

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

        // experimental!!!
        // Equality constraint (w1 + w2 + w3 = 1)
        _E.resize(25, 1);
        _E.setConstant(ukfZero);
        _e.resize(1);

        _E(21, 0) = _E(22, 0) = _E(23, 0) = -1.0;
        _e(0) = 1;
    }

    virtual ~Ridg_BiExp_FW()
    {
    }

    virtual void F(ukfMatrixType & /* X */, ukfVectorType /* s */, const ukfMatrixType & /* &covMatrix */) const;
    ukfPrecisionType cosine_similarity(vec3_t &First, vec3_t &Second) const;
    virtual void F(ukfMatrixType & /* X */) const;
    virtual void H(const ukfMatrixType &X, ukfMatrixType &Y) const;

    virtual void State2Tensor3T(const State &x, const vec3_t &old_m, vec3_t &m1, vec3_t &l1, vec3_t &m2, vec3_t &l2, vec3_t &m3, vec3_t &l3);
    virtual void State2Tensor3T(const State &x, const vec3_t &old_m, vec3_t &m1, vec3_t &m2, vec3_t &m3);

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
    ukfMatrixType &Q;
    ukfMatrixType &fcs;
    ukfMatrixType &nu;
    vector<vector<unsigned>> &conn;

    signalMaskType &signal_mask;

    const ukfPrecisionType fista_lambda;
    const ukfPrecisionType max_odf_thresh;
};

#endif //RIDG_BiExp_FW__