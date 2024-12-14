/**
 * @file BInterpolate.tpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief This is the implementation of BInterpolate
 * @version 0.1
 * @date 2024-11-24
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef BINTERPOLATE_TPP
#define BINTERPOLATE_TPP

#include "BInterpolate.hpp"

template <int N, typename Real>
BInterpolate<N, Real>::BInterpolate(const std::vector<Real> &t,
                                    const std::vector<Real> &y,
                                    const int               &method,
                                    const std::vector<Real> &boundary_condition,
                                    const int                check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<Real>();
    if (check == 1)
    {
        if (!std::is_sorted(_t.begin(), _t.end()))
        {
            int              t_size = _t.size();
            std::vector<int> idx(t_size);

            for (int i = 0; i < t_size; ++i)
            {
                idx[i] = i;
            }
            std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });

            for (int i = 0; i < t_size; ++i)
            {
                _t[i] = t[idx[i]];
                _y[i] = y[idx[i]];
            }
        }
    }
    interpolate(t, y, method, boundary_condition);
}

template <int N, typename Real>
BInterpolate<N, Real>::BInterpolate(const Json::Value &json)
{
    _bspline = BSpline<Real>();
    // load from json file
    _method   = json.isMember("method") ? json["method"].asInt() : 0;
    int check = json.isMember("check") ? json["check"].asInt() : 0;

    if (json.isMember("t"))
    {
        const Json::Value &t_json = json["t"];
        for (auto &t : t_json)
        {
            _t.push_back(Real(t.asDouble()));
        }
    }
    if (json.isMember("y"))
    {
        const Json::Value &y_json = json["y"];
        for (auto &y : y_json)
        {
            _y.push_back(Real(y.asDouble()));
        }
    }
    if (json.isMember("boundary_condition"))
    {
        const Json::Value &boundary_condition_json = json["boundary_condition"];
        for (auto &boundary_condition : boundary_condition_json)
        {
            _boundary_condition.push_back(Real(boundary_condition.asDouble()));
        }
    }
    else
    {
        _boundary_condition = std::vector<Real>(N, 0.0);
    }

    if (check == 1)
    {
        if (!std::is_sorted(_t.begin(), _t.end()))
        {
            std::vector<Real> t      = _t;
            std::vector<Real> y      = _y;
            int               t_size = t.size();
            std::vector<int>  idx(t_size);

            for (int i = 0; i < t_size; ++i)
            {
                idx[i] = i;
            }
            std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });

            for (int i = 0; i < t_size; ++i)
            {
                _t[i] = t[idx[i]];
                _y[i] = y[idx[i]];
            }
        }
    }
    interpolate(_t, _y, _method, _boundary_condition);
}

template <int N, typename Real>
void
BInterpolate<N, Real>::interpolate(const std::vector<Real> &t,
                                   const std::vector<Real> &y,
                                   const int               &method, // 0 for periodic, 1 for complete,
                                                                    // 2 for natural, 3 for not-a-knot.
                                   const std::vector<Real> &boundary_condition)
{
    /**
     * @details Here are the details for the interpolation.
     *
     * For B-Spline, the interpolation is easy, we can easily construct the matrix A
     * and vector b.
     *
     * Firstly, for the \f$f(x_i) = y_i\f$, we just calculate the basis function that is non-zero at x_i,
     * and let them be the coefficients.
     *
     */
    if (method == 0)
    {
        /**
         * @details **Periodic Condition**
         *
         * For the periodic condition, we can have \f$n-1\f$ more condtions that is
         * \f[
         * \frac{\mathrm{d}^i f(x_1)}{\mathrm{d} x^i} = \frac{\mathrm{d}^i f(x_n)}{\mathrm{d} x^i}
         * \f]
         *
         * Recall the form of B-Spline, then we can easily construct the matrix A and vector b.
         *
         */
        int                                                 t_size = t.size();
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A(t_size + N - 1, t_size + N - 1);
        Eigen::Matrix<Real, Eigen::Dynamic, 1>              b(t_size + N - 1);
        std::vector<Real>                                   tmp_coeff(t_size + N - 1, 1.0);
        BSpline<Real>                                       tmp_spline(tmp_coeff, t, N);
        for (int i = 0; i < N - 1; i++)
        {
            std::vector<Real> diff_basis1 = tmp_spline.basis_derivative(t[0], i + 1);
            std::vector<Real> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], i + 1);
            for (int j = 0; j < N; ++j)
            {
                A(i, j)              = diff_basis1[j];
                A(i, t_size + j - 1) = -diff_basis2[j + 1];
            }
            b(i) = Real(0.0);
        }

        for (int i = 0; i < t_size; ++i)
        {
            std::vector<Real> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                for (int j = 0; j < N; ++j)
                {
                    A(N - 1 + i, i + j) = basis[j];
                }
            }
            else
            {
                for (int j = 0; j < N; ++j)
                {
                    A(N - 1 + i, i + j) = basis[j + 1];
                }
            }
            b(N - 1 + i) = y[i];
        }

        Eigen::PartialPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> solver;
        solver.compute(A);

        Eigen::Matrix<Real, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<Real>                      coeffs(x.data(), x.data() + x.size());

        _bspline = BSpline<Real>(coeffs, t, N);
        return;
    }
    else if (method == 1)
    {
        /**
         * @details **Complete Condition**
         *
         * This condition is based on all the derivates at the beginning point are given.
         *
         * The more \f$n-1\f$ equations is easy to construct.
         *
         */
        int t_size = t.size();

        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A_first(N, N);
        Eigen::Matrix<Real, Eigen::Dynamic, 1>              b_first(N);
        std::vector<Real>                                   tmp_coeff(t_size + N - 1, 1.0);
        BSpline<Real>                                       tmp_spline(tmp_coeff, t, N);
        for (int i = 0; i < N - 1; i++)
        {
            std::vector<Real> diff_basis = tmp_spline.basis_derivative(t[0], i + 1);
            for (int j = 0; j < N; ++j)
            {
                A_first(i, j) = diff_basis[j];
            }
            b_first(i) = boundary_condition[i];
        }

        Eigen::PartialPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> solver;
        solver.compute(A_first);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> x_first = solver.solve(b_first);

        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A(t_size - 1, t_size - 1);
        Eigen::Matrix<Real, Eigen::Dynamic, 1>              b(t_size - 1);
        for (int i = 0; i < N - 1 && i < t_size - 1; i++)
        {
            std::vector<Real> basis = tmp_spline.get_basis(t[i + 1]);
            for (int j = 0; j <= i; j++)
            {
                A(i, j) = basis[N - i + j];
            }
            b(i) = y[i + 1];
            for (int j = 0; j < N - i - 1; j++)
            {
                b(i) -= basis[j + 1] * x_first(i + j);
            }
        }
        for (int i = N - 1; i < t_size - 1; i++)
        {
            std::vector<Real> basis = tmp_spline.get_basis(t[i + 1]);
            for (int j = 0; j < N; j++)
            {
                A(i, j) = basis[j];
            }
            b(i) = y[i + 1];
        }

        solver.compute(A);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> x = solver.solve(b);
        // add first to the begin of x
        Eigen::Matrix<Real, Eigen::Dynamic, 1> x_new(N + t_size - 1);
        x_new << x_first, x;
        std::vector<Real> coeff(x_new.data(), x_new.data() + x_new.size());
        _bspline = BSpline<Real>(coeff, t, N);
        return;

        // Eigen::SparseMatrix<Real, Eigen::RowMajor> A(t_size + N - 1, t_size +
        // N - 1); Eigen::MatrixXd A(t_size+N-1, t_size+N-1); Eigen::VectorXd
        // b(t_size+N-1);
        // Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A(t_size + N - 1, t_size + N - 1);
        // Eigen::Matrix<Real, Eigen::Dynamic, 1>              b(t_size + N - 1);
        // // std::vector<Eigen::Triplet<Real>>      triplets;
        // // triplets.reserve(N * (t_size + N - 1));
        // std::vector<Real> tmp_coeff(t_size + N - 1, 1.0);
        // BSpline<Real>     tmp_spline(tmp_coeff, t, N);

        // for (int i = 0; i < N - 1; i++)
        // {
        //     std::vector<Real> diff_basis = tmp_spline.basis_derivative(t[0], i + 1);
        //     for (int j = 0; j < N; ++j)
        //     {
        //         // triplets.push_back(Eigen::Triplet<Real>(i, j,
        //         // diff_basis[j]));
        //         A(i, j) = diff_basis[j];
        //     }
        //     b(i) = boundary_condition[i];
        // }

        // for (int i = 0; i < t_size; ++i)
        // {
        //     std::vector<Real> basis = tmp_spline.get_basis(t[i]);
        //     if (i == 0)
        //     {
        //         for (int j = 0; j < N; ++j)
        //         {
        //             // triplets.push_back(Eigen::Triplet<Real>(N - 1 + i, i + j,
        //             // basis[j]));
        //             A(N - 1 + i, i + j) = basis[j];
        //         }
        //     }
        //     else
        //     {
        //         for (int j = 0; j < N; ++j)
        //         {
        //             // triplets.push_back(Eigen::Triplet<Real>(N - 1 + i, i + j,
        //             // basis[j + 1]));
        //             A(N - 1 + i, i + j) = basis[j + 1];
        //         }
        //     }
        //     b(N - 1 + i) = y[i];
        // }

        // // A.setFromTriplets(triplets.begin(), triplets.end());
        // // A.makeCompressed();
        // // 计算矩阵条件数
        // // use more stable solver, not LU, iteration method
        // // use the most accurate solver
        // // Eigen::PartialPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> solver;
        // Eigen::FullPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> solver;
        // solver.compute(A);
        // // if (solver.info() != Eigen::Success)
        // // {
        // //     std::cerr << "Failed to factorize the matrix" << std::endl;
        // // }
        // // Eigen::VectorXd x = solver.solve(b);
        // Eigen::Matrix<Real, Eigen::Dynamic, 1> x = solver.solve(b);

        // // if (solver.info() != Eigen::Success)
        // // {
        // //     std::cerr << "Failed to solve the system of equations" <<
        // //     std::endl;
        // // }

        // std::vector<Real> coeffs(x.data(), x.data() + x.size());

        // // std::cout << "A * x - b = " << std::endl << A * x - b << std::endl;
        // _bspline = BSpline<Real>(coeffs, t, N);
        // return;
    }
    else
    {
        throw std::invalid_argument("method must be 0,1");
        return;
    }
}

template <int N, typename Real>
Real
BInterpolate<N, Real>::operator()(const Real x)
{
    return _bspline(x);
}

template <int N, typename Real>
BSpline<Real>
BInterpolate<N, Real>::getBSpline() const
{
    return _bspline;
}

template <int N, typename Real>
Real
BInterpolate<N, Real>::derivative(const Real x, const int n)
{
    return _bspline.derivative(x, n);
}

#endif // BINTERPOLATE_TPP