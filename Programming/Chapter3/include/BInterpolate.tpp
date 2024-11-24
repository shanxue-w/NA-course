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
            std::sort(idx.begin(),
                      idx.end(),
                      [&](int i, int j) { return t[i] < t[j]; });

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
void
BInterpolate<N, Real>::interpolate(
    const std::vector<Real> &t,
    const std::vector<Real> &y,
    const int               &method, // 0 for periodic, 1 for complete,
                                     // 2 for natural, 3 for not-a-knot.
    const std::vector<Real> &boundary_condition)
{
    // (void)method;
    if (method == 0)
    {
        int                                                 t_size = t.size();
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A(t_size + N - 1,
                                                              t_size + N - 1);
        Eigen::Matrix<Real, Eigen::Dynamic, 1>              b(t_size + N - 1);
        std::vector<Real> tmp_coeff(t_size + N - 1, 1.0);
        BSpline<Real>     tmp_spline(tmp_coeff, t, N);
        for (int i = 0; i < N - 1; i++)
        {
            std::vector<Real> diff_basis1 =
                tmp_spline.basis_derivative(t[0], i + 1);
            std::vector<Real> diff_basis2 =
                tmp_spline.basis_derivative(t[t_size - 1], i + 1);
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

        Eigen::PartialPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
            solver;
        solver.compute(A);

        Eigen::Matrix<Real, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<Real> coeffs(x.data(), x.data() + x.size());

        _bspline = BSpline<Real>(coeffs, t, N);
        return;
    }
    else if (method == 1)
    {
        int t_size = t.size();
        // Eigen::SparseMatrix<Real, Eigen::RowMajor> A(t_size + N - 1, t_size +
        // N - 1); Eigen::MatrixXd A(t_size+N-1, t_size+N-1); Eigen::VectorXd
        // b(t_size+N-1);
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A(t_size + N - 1,
                                                              t_size + N - 1);
        Eigen::Matrix<Real, Eigen::Dynamic, 1>              b(t_size + N - 1);
        // std::vector<Eigen::Triplet<Real>>      triplets;
        // triplets.reserve(N * (t_size + N - 1));
        std::vector<Real> tmp_coeff(t_size + N - 1, 1.0);
        BSpline<Real>     tmp_spline(tmp_coeff, t, N);

        for (int i = 0; i < N - 1; i++)
        {
            std::vector<Real> diff_basis =
                tmp_spline.basis_derivative(t[0], i + 1);
            for (int j = 0; j < N; ++j)
            {
                // triplets.push_back(Eigen::Triplet<Real>(i, j,
                // diff_basis[j]));
                A(i, j) = diff_basis[j];
            }
            b(i) = boundary_condition[i];
        }

        for (int i = 0; i < t_size; ++i)
        {
            std::vector<Real> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                for (int j = 0; j < N; ++j)
                {
                    // triplets.push_back(Eigen::Triplet<Real>(N - 1 + i, i + j,
                    // basis[j]));
                    A(N - 1 + i, i + j) = basis[j];
                }
            }
            else
            {
                for (int j = 0; j < N; ++j)
                {
                    // triplets.push_back(Eigen::Triplet<Real>(N - 1 + i, i + j,
                    // basis[j + 1]));
                    A(N - 1 + i, i + j) = basis[j + 1];
                }
            }
            b(N - 1 + i) = y[i];
        }

        // A.setFromTriplets(triplets.begin(), triplets.end());
        // A.makeCompressed();
        // 计算矩阵条件数
        // use more stable solver, not LU, iteration method
        // use the most accurate solver
        Eigen::PartialPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
            solver;
        solver.compute(A);
        // if (solver.info() != Eigen::Success)
        // {
        //     std::cerr << "Failed to factorize the matrix" << std::endl;
        // }
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> x = solver.solve(b);

        // if (solver.info() != Eigen::Success)
        // {
        //     std::cerr << "Failed to solve the system of equations" <<
        //     std::endl;
        // }

        std::vector<Real> coeffs(x.data(), x.data() + x.size());

        // std::cout << "A * x - b = " << std::endl << A * x - b << std::endl;
        _bspline = BSpline<Real>(coeffs, t, N);
        return;
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