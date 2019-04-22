#include "filter_ridg.h"

// TODO: Check with model on my notebook and adjust model in H function

// 2T Bi-Exponential model with spherical ridgelets //
// Functions for 2-tensor bi-exponential simple model.
void Ridg_BiExp_FW::F(ukfMatrixType &X) const
{
    assert(_signal_dim > 0);
    assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
           (X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
            X.cols() == 1));

    // Need to add ridgelets here

    for (unsigned int i = 0; i < X.cols(); ++i)
    {
        // Normalize the direction vectors.

        // Tensor 1
        ukfPrecisionType norm_inv = ukfZero; // 1e-16;
        norm_inv += X(0, i) * X(0, i);
        norm_inv += X(1, i) * X(1, i);
        norm_inv += X(2, i) * X(2, i);

        norm_inv = ukfOne / sqrt(norm_inv);
        X(0, i) *= norm_inv;
        X(1, i) *= norm_inv;
        X(2, i) *= norm_inv;

        // Tensor 2
        norm_inv = ukfZero; // 1e-16;
        norm_inv += X(7, i) * X(7, i);
        norm_inv += X(8, i) * X(8, i);
        norm_inv += X(9, i) * X(9, i);

        norm_inv = ukfOne / sqrt(norm_inv);
        X(7, i) *= norm_inv;
        X(8, i) *= norm_inv;
        X(9, i) *= norm_inv;

        // Tensor 3
        norm_inv = ukfZero; // 1e-16;
        norm_inv += X(14, i) * X(14, i);
        norm_inv += X(15, i) * X(15, i);
        norm_inv += X(16, i) * X(16, i);

        norm_inv = ukfOne / sqrt(norm_inv);
        X(14, i) *= norm_inv;
        X(15, i) *= norm_inv;
        X(16, i) *= norm_inv;

        // Check that the eigenvalues are greater or equal to the minimum value
        // and less or equal to the maximum value
        // Tensor 1
        X(3, i) = std::max(X(3, i), _lambda_min_fast_diffusion);
        X(4, i) = std::max(X(4, i), _lambda_min_fast_diffusion);
        X(5, i) = std::max(X(5, i), _lambda_min_slow_diffusion);
        X(6, i) = std::max(X(6, i), _lambda_min_slow_diffusion);

        X(3, i) = std::min(X(3, i), _lambda_max_diffusion);
        X(4, i) = std::min(X(4, i), _lambda_max_diffusion);
        X(5, i) = std::min(X(5, i), _lambda_max_diffusion);
        X(6, i) = std::min(X(6, i), _lambda_max_diffusion);

        // Tensor 2
        X(10, i) = std::max(X(10, i), _lambda_min_fast_diffusion);
        X(11, i) = std::max(X(11, i), _lambda_min_fast_diffusion);
        X(12, i) = std::max(X(12, i), _lambda_min_slow_diffusion);
        X(13, i) = std::max(X(13, i), _lambda_min_slow_diffusion);

        X(10, i) = std::min(X(10, i), _lambda_max_diffusion);
        X(11, i) = std::min(X(11, i), _lambda_max_diffusion);
        X(12, i) = std::min(X(12, i), _lambda_max_diffusion);
        X(13, i) = std::min(X(13, i), _lambda_max_diffusion);

        // Tensor 3
        X(17, i) = std::max(X(17, i), _lambda_min_fast_diffusion);
        X(18, i) = std::max(X(18, i), _lambda_min_fast_diffusion);
        X(19, i) = std::max(X(19, i), _lambda_min_slow_diffusion);
        X(20, i) = std::max(X(20, i), _lambda_min_slow_diffusion);

        X(10, i) = std::min(X(17, i), _lambda_max_diffusion);
        X(11, i) = std::min(X(18, i), _lambda_max_diffusion);
        X(12, i) = std::min(X(19, i), _lambda_max_diffusion);
        X(13, i) = std::min(X(20, i), _lambda_max_diffusion);

        // Free water
        X(23, i) = CheckZero(X(14, i));
    }
}

void Ridg_BiExp_FW::H(const ukfMatrixType &X,
                      ukfMatrixType &Y) const
{
    assert(_signal_dim > 0);
    assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
           (X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
            X.cols() == 1));
    assert(Y.rows() == static_cast<unsigned int>(_signal_dim) &&
           (Y.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
            Y.cols() == 1));
    assert(_signal_data);

    const stdVec_t &gradients = _signal_data->gradients();
    const ukfVectorType &b = _signal_data->GetBValues();
    for (unsigned int i = 0; i < X.cols(); ++i)
    {
        // Normalize directions.
        vec3_t m1;
        initNormalized(m1, X(0, i), X(1, i), X(2, i));
        vec3_t m2;
        initNormalized(m2, X(7, i), X(8, i), X(9, i));
        vec3_t m3;
        initNormalized(m3, X(14, i), X(15, i), X(16, i));

        // Tensor 1 lambdas
        ukfPrecisionType l11 = std::max(X(3, i), _lambda_min_fast_diffusion);
        ukfPrecisionType l12 = std::max(X(4, i), _lambda_min_fast_diffusion);
        ukfPrecisionType l13 = std::max(X(5, i), _lambda_min_slow_diffusion);
        ukfPrecisionType l14 = std::max(X(6, i), _lambda_min_slow_diffusion);

        l11 = std::min(l11, _lambda_max_diffusion);
        l12 = std::min(l12, _lambda_max_diffusion);
        l13 = std::min(l13, _lambda_max_diffusion);
        l14 = std::min(l14, _lambda_max_diffusion);

        // Tensor 2 lambdas
        ukfPrecisionType l21 = std::max(X(10, i), _lambda_min_fast_diffusion);
        ukfPrecisionType l22 = std::max(X(11, i), _lambda_min_fast_diffusion);
        ukfPrecisionType l23 = std::max(X(12, i), _lambda_min_slow_diffusion);
        ukfPrecisionType l24 = std::max(X(13, i), _lambda_min_slow_diffusion);

        l21 = std::min(l21, _lambda_max_diffusion);
        l22 = std::min(l22, _lambda_max_diffusion);
        l23 = std::min(l23, _lambda_max_diffusion);
        l24 = std::min(l24, _lambda_max_diffusion);

        // Tensor 3 lambdas
        ukfPrecisionType l31 = std::max(X(17, i), _lambda_min_fast_diffusion);
        ukfPrecisionType l32 = std::max(X(18, i), _lambda_min_fast_diffusion);
        ukfPrecisionType l33 = std::max(X(19, i), _lambda_min_slow_diffusion);
        ukfPrecisionType l34 = std::max(X(20, i), _lambda_min_slow_diffusion);

        l31 = std::min(l31, _lambda_max_diffusion);
        l32 = std::min(l32, _lambda_max_diffusion);
        l33 = std::min(l33, _lambda_max_diffusion);
        l34 = std::min(l34, _lambda_max_diffusion);

        // Flip if necessary.
        if (m1[0] < 0)
        {
            m1 = -m1;
        }
        if (m2[0] < 0)
        {
            m2 = -m2;
        }
        if (m3[0] < 0)
        {
            m3 = -m3;
        }

        // Get free water weight from state
        const ukfPrecisionType w = CheckZero(X(23, i));

        // Fill in lambdas matricies
        diagmat3_t lambdas11, lambdas12, lambdas21, lambdas22, lambdas31, lambdas32;
        lambdas11.diagonal()[0] = l11;
        lambdas11.diagonal()[1] = l12;
        lambdas11.diagonal()[2] = l12;

        lambdas12.diagonal()[0] = l13;
        lambdas12.diagonal()[1] = l14;
        lambdas12.diagonal()[2] = l14;

        lambdas21.diagonal()[0] = l21;
        lambdas21.diagonal()[1] = l22;
        lambdas21.diagonal()[2] = l22;

        lambdas22.diagonal()[0] = l23;
        lambdas22.diagonal()[1] = l24;
        lambdas22.diagonal()[2] = l24;

        lambdas31.diagonal()[0] = l31;
        lambdas31.diagonal()[1] = l32;
        lambdas31.diagonal()[2] = l32;

        lambdas32.diagonal()[0] = l33;
        lambdas32.diagonal()[1] = l34;
        lambdas32.diagonal()[2] = l34;

        // Calculate diffusion matrix.
        const mat33_t &D1 = diffusion(m1, lambdas11);
        const mat33_t &D1t = diffusion(m1, lambdas12);
        const mat33_t &D2 = diffusion(m2, lambdas21);
        const mat33_t &D2t = diffusion(m2, lambdas22);
        const mat33_t &D3 = diffusion(m3, lambdas31);
        const mat33_t &D3t = diffusion(m3, lambdas32);

        // Reconstruct signal by the means of the model
        for (int j = 0; j < _signal_dim; ++j)
        {
            // u = gradient direction considered
            const vec3_t &u = gradients[j];

            Y(j, i) =
                (1 - w) * (weights_on_tensors_[0] * (_w_fast_diffusion * std::exp(-b[j] * u.dot(D1 * u)) + (1 - _w_fast_diffusion) * std::exp(-b[j] * u.dot(D1t * u))) +
                 weights_on_tensors_[1] * (_w_fast_diffusion * std::exp(-b[j] * u.dot(D2 * u)) + (1 - _w_fast_diffusion) * std::exp(-b[j] * u.dot(D2t * u))) + 
                 weights_on_tensors_[2] * (_w_fast_diffusion * std::exp(-b[j] * u.dot(D3 * u)) + (1 - _w_fast_diffusion) * std::exp(-b[j] * u.dot(D3t * u)))) + 
                 w * std::exp(-b[j] * u.dot(m_D_iso * u));
        }
    }
}

void Ridg_BiExp_FW::State2Tensor3T(const State &x, const vec3_t &old_m, vec3_t &m1, vec3_t &l1, vec3_t &m2, vec3_t &l2, vec3_t &m3, vec3_t &l3)
{
    // Orientations;
    initNormalized(m1, x[0], x[1], x[2]);
    initNormalized(m2, x[7], x[8], x[9]);
    initNormalized(m3, x[14], x[15], x[16]);

    // Tensor 1
    // Lambda fast diffusion
    l1[0] = std::max(x(3), _lambda_min_fast_diffusion);
    l1[1] = std::max(x(4), _lambda_min_fast_diffusion);
    l1[2] = l1[1];

    // Tensor 2
    // Lambda fast diffusion
    l2[0] = std::max(x(10), _lambda_min_fast_diffusion);
    l2[1] = std::max(x(11), _lambda_min_fast_diffusion);
    l2[2] = l2[1];

    // Tensor 3
    // Lambda fast diffusion
    l3[0] = std::max(x(17), _lambda_min_fast_diffusion);
    l3[1] = std::max(x(18), _lambda_min_fast_diffusion);
    l3[2] = l3[1];

    // Flip orientations if necessary. (For m1 it should not happen, maybe for
    // m2.)
    if (m1[0] * old_m[0] + m1[1] * old_m[1] + m1[2] * old_m[2] < 0)
    {
        m1 = -m1;
    }
    if (m2[0] * old_m[0] + m2[1] * old_m[1] + m2[2] * old_m[2] < 0)
    {
        m2 = -m2;
    }
    if (m3[0] * old_m[0] + m3[1] * old_m[1] + m3[2] * old_m[2] < 0)
    {
        m3 = -m3;
    }
}
