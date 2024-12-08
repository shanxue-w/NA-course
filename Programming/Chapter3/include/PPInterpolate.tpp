/**
 * @file PPInterpolate.tpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implmentation of PP-form interpolation
 * @version 0.1
 * @date 2024-11-24
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "PPInterpolate.hpp"
#include <algorithm>
#include <iostream>

template <int N, typename Real>
PPInterpolate<N, Real>::PPInterpolate(const std::vector<Real> &t, // nodes
                                      const std::vector<Real> &y, // values
                                      const int method,           // 0 for periodic, 1 for complete,
                                      // 2 for natural, 3 for not-a-knot.
                                      const std::vector<Real> &boundary_condition,
                                      const int check) // whether to check the order of t
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    poly = PPoly<Real>();
    if (check)
    {
        if (!std::is_sorted(t.begin(), t.end()))
        {
            std::vector<int> idx(t.size());

            for (size_t i = 0; i < t.size(); ++i)
            {
                idx[i] = i;
            }
            std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });
            // std::vector<Real> t_sorted(t.size()), y_sorted(t.size());

            for (size_t i = 0; i < t.size(); ++i)
            {
                _t[i] = t[idx[i]];
                _y[i] = y[idx[i]];
            }
        }
    }
    interpolate(t, y, method, boundary_condition); // interpolate
}

template <int N, typename Real>
void
PPInterpolate<N, Real>::interpolate(
    const std::vector<Real> &t,      // nodes
    const std::vector<Real> &y,      // values
    const int                method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                     // not-a-knot.
    const std::vector<Real> &boundary_condition) // boundary condition
{
    /**
     * @details Here are the details of the interpolation:
     * Let the polynomial for the interval [x_i, x_{i+1}] be defined as:
     * \f[
     * p_i(x) = y_i + C_{i,1}(x-x_i) + \cdots + C_{i,N}(x-x_i)^N
     * \f]
     * and for the next interval \f$[x_{i+1}, x_{i+2}]\f$ as:
     * \f[
     * p_{i+1}(x) = y_{i+1} + C_{i+1,1}(x-x_{i+1}) + \cdots +
     * C_{i+1,N}(x-x_{i+1})^N
     * \f]
     *
     * By requiring continuity at the boundary \f$x_{i+1}\f$, i.e.,
     * \f$p_i(x_{i+1}) = p_{i+1}(x_{i+1})\f$, define: \f[ \Delta x_i = x_{i+1} -
     * x_i, \quad \Delta y_i = y_{i+1} - y_i \f]
     *
     * The coefficient \f$C_{i,N}\f$ is then computed as:
     * \f[
     * C_{i,N} = \frac{\Delta y_i - C_{i,1}\Delta x_i - \cdots - C_{i,N-1}\Delta
     * x_i^{N-1}}{\Delta x_i^N} \f]
     *
     * By matching the \f$j\f$-th derivatives at \f$x_{i+1}\f$, i.e.,
     * \f$p_i^{(j)}(x_{i+1}) = p_{i+1}^{(j)}(x_{i+1})\f$, we derive:
     * \f[
     * C_{i+1,j} = \sum_{k=0}^{N-j} C_{j+k}^{k} C_{i,j+k} (\Delta x_i)^k
     * \f]
     *
     * Replacing \f$C_{i,N}\f$ in the above equation, we get:
     * \f[
     * C_{i+1,j} = \sum_{k=0}^{N-j-1} C_{j+k}^{k} C_{i,j+k} (\Delta x_i)^k
     * + C_{N}^{N-j} \frac{\Delta y_i - C_{i,1}\Delta x_i - \cdots -
     * C_{i,N-1}\Delta x_i^{N-1}}{\Delta x_i^j}
     * \f]
     *
     * The matrix equation is then:
     * \f$ b + A C_i = C_{i+1} \f$
     */
    Eigen::Matrix<Real, N, N> CMatrix = Generate_Cmatrix();
    // std::cout << "CMatrix: " << CMatrix << std::endl;
    if (method == 0)
    {
        int                            t_size = t.size();
        std::vector<std::vector<Real>> coeffs(
            t_size - 1, std::vector<Real>(N)); // coefficients of the polynomial
        Eigen::Matrix<Real, N - 1, 1>     C = Eigen::Matrix<Real, N - 1, 1>::Zero(N - 1);
        Eigen::Matrix<Real, N - 1, N - 1> A_ahead =
            Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1);
        Eigen::Matrix<Real, N - 1, 1>     b_ahead = Eigen::Matrix<Real, N - 1, 1>::Zero(N - 1);
        Eigen::Matrix<Real, N - 1, N - 1> A_behind =
            Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1);
        Eigen::Matrix<Real, N - 1, 1> b_behind  = Eigen::Matrix<Real, N - 1, 1>::Zero(N - 1);
        Eigen::Matrix<Real, N, 1>     tmp_lists = Eigen::Matrix<Real, N, 1>::Ones(N);
        int                           ahead = 0, behind = t_size - 1;

        Eigen::Matrix<Real, N - 1, N - 1> A_copy =
            Eigen::Matrix<Real, N - 1, N - 1>::Zero(N - 1, N - 1);
        Eigen::Matrix<Real, N - 1, 1> b_copy = Eigen::Matrix<Real, N - 1, 1>::Zero(N - 1);
        Real                          dx     = t[ahead + 1] - t[ahead];
        Real                          dy     = y[ahead + 1] - y[ahead];
        Real                          tmpx   = dx;
        Real                          C_iN   = dy; // C_{i,N}
        std::cout << "CMatrix: " << CMatrix << std::endl;
        while (ahead < behind)
        {
            Real max_value1 = A_ahead.cwiseAbs().maxCoeff();
            Real max_value2 = A_behind.cwiseAbs().maxCoeff();
            if (max_value1 < max_value2)
            {
                // first compute the ahead step
                // C_{i+1} = A_i C_i + b_i
                // C_i = A_ahead_i C_0 + b_ahead_i
                // A_ahead_{i+1} = A_i A_ahead_i
                // b_ahead_{i+1} = A_i b_ahead_i + b_i
                tmp_lists    = Eigen::Matrix<Real, N, 1>::Ones(N);
                tmp_lists(0) = dx;
                for (int j = 1; j < N; j++)
                {
                    tmp_lists(j) = tmp_lists(j - 1) * dx;
                }
                // compute forward A_ahead, b_ahead
                A_copy = Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1);
                for (int j = 0; j < N - 1; j++)
                {
                    for (int k = j + 1; k < N - 1; k++)
                    {
                        A_copy(j, k) = CMatrix(j, k) * tmp_lists(k - j - 1);
                    }
                    Real tmp  = CMatrix(j, N - 1) / tmp_lists(j);
                    b_copy(j) = tmp * dy;
                    for (int k = 0; k < N - 1; k++)
                    {
                        A_copy(j, k) -= tmp * tmp_lists(k);
                    }
                }

                A_ahead = A_copy * A_ahead;
                b_ahead = A_copy * b_ahead + b_copy;
                ahead++;
            }
            else
            {
                // then compute the backward step
                // C_{i+1} = A C_i + b
                // C_i = A^{-1} C_{i+1} - A^{-1} b
                // C_{i+1} = A_behind C_N + b_behind
                // A_behind_new = A^{-1} A_behind
                // b_behind_new = A^{-1} b_behind - A^{-1} b
                dx           = t[behind] - t[behind - 1];
                dy           = y[behind] - y[behind - 1];
                tmpx         = dx;
                tmp_lists    = Eigen::Matrix<Real, N, 1>::Ones(N);
                tmp_lists(0) = dx;
                for (int j = 1; j < N; j++)
                {
                    tmp_lists(j) = tmp_lists(j - 1) * dx;
                }
                // compute backward A_behind, b_behind
                A_copy = Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1);
                for (int j = 0; j < N - 1; j++)
                {
                    for (int k = j + 1; k < N - 1; k++)
                    {
                        A_copy(j, k) = CMatrix(j, k) * tmp_lists(k - j - 1);
                    }
                    Real tmp  = CMatrix(j, N - 1) / tmp_lists(j);
                    b_copy(j) = tmp * dy;
                    for (int k = 0; k < N - 1; k++)
                    {
                        A_copy(j, k) -= tmp * tmp_lists(k);
                    }
                }
                // A_copy   = A_copy.inverse();
                // Eigen::Matrix<Real, N - 1, N - 1> A_copy_inverse = A_copy.inverse();
                // A_behind                                         = A_copy_inverse * A_behind;
                // b_behind = A_copy_inverse * (b_behind - b_copy);
                // 通过ColPivHouseholderQR求解线性系统 A_copy * x = b_behind
                Eigen::ColPivHouseholderQR<Eigen::Matrix<Real, N - 1, N - 1>> solver;
                solver.compute(A_copy);            // 计算QR分解
                A_behind = solver.solve(A_behind); // 解线性方程 A_copy * x = A_behind

                // 求解b_behind
                b_behind =
                    solver.solve(b_behind - b_copy); // 解线性方程 A_copy * x = b_behind - b_copy
                behind--;
            }

            // find the maximum
            max_value1 = A_ahead.cwiseAbs().maxCoeff();
            max_value2 = A_behind.cwiseAbs().maxCoeff();
            // if (max_value1 > 1e6 and max_value2 > 1e6 and max_value3 > 1e6 and max_value4 > 1e6)
            // {
            //     A_ahead /= Real(1e4);
            //     b_ahead /= Real(1e4);
            //     A_behind /= Real(1e4);
            //     b_behind /= Real(1e4);
            // }

            if (ahead == behind)
            {
                break;
            }
        }
        Eigen::Matrix<Real, N - 1, N - 1> A = A_ahead - A_behind;
        Eigen::Matrix<Real, N - 1, 1>     b = -(b_ahead - b_behind);
        // Real                              maxvalue = A.cwiseAbs().maxCoeff() / Real(100);
        // A /= maxvalue;
        // b /= maxvalue;
        // Eigen::ColPivHouseholderQR<Eigen::Matrix<Real, N - 1, N - 1>> solver;
        Eigen::FullPivHouseholderQR<Eigen::Matrix<Real, N - 1, N - 1>> solver;
        // Eigen::PartialPivLU<Eigen::Matrix<Real, N - 1, N - 1>> solver;
        solver.compute(A);
        C                                 = solver.solve(b);
        Eigen::Matrix<Real, N - 1, 1> C_0 = C; // save the initial value of C

        dx   = t[1] - t[0];
        dy   = y[1] - y[0];
        tmpx = dx;
        C_iN = dy;
        for (int i = 0; i < N - 1; i++)
        {
            C_iN -= C(i) * tmpx;
            tmpx *= dx;
        }
        C_iN /= tmpx;
        std::vector<Real> tmp_coeffs(N, 0.0);
        tmp_coeffs[0] = y[0];
        for (int i = 0; i < N - 1; i++)
        {
            tmp_coeffs[i + 1] = C(i);
        }
        tmp_coeffs.push_back(C_iN);
        coeffs[0] = tmp_coeffs;

        for (int i = 0; i < ahead / 2; i++)
        {
            std::cout << "i = " << i << std::endl;
            std::cout << "C = " << C << std::endl;
            dx           = t[i + 1] - t[i];
            dy           = y[i + 1] - y[i];
            tmp_lists(0) = dx;
            for (int j = 1; j < N; j++)
            {
                tmp_lists(j) = tmp_lists(j - 1) * dx;
            }

            A = Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1);
            for (int j = 0; j < N - 1; j++)
            {
                for (int k = j + 1; k < N - 1; k++)
                {
                    A(j, k) = CMatrix(j, k) * tmp_lists(k - j - 1);
                }
                Real tmp = CMatrix(j, N - 1) / tmp_lists(j);
                b(j)     = tmp * dy;
                for (int k = 0; k < N - 1; k++)
                {
                    A(j, k) -= tmp * tmp_lists(k);
                }
            }
            C = b + A * C;
            // Eigen::Matrix<Real, N-1, 1> D = b + A * C;
            // use QR to solve D

            dx   = t[i + 2] - t[i + 1];
            dy   = y[i + 2] - y[i + 1];
            tmpx = dx;
            C_iN = dy;
            for (int k = 0; k < N - 1; k++)
            {
                C_iN -= C(k) * tmpx;
                tmpx *= dx;
            }
            C_iN /= tmpx;
            tmp_coeffs[0] = y[i + 1];
            for (int k = 0; k < N - 1; k++)
            {
                tmp_coeffs[k + 1] = C(k);
            }
            tmp_coeffs[N] = C_iN;
            coeffs[i + 1] = tmp_coeffs;
        }
        C = C_0;
        for (int i = t_size - 1; i > ahead / 2; i--)
        {
            std::cout << "i = " << i << std::endl;
            std::cout << "C = " << C << std::endl;
            dx           = t[i] - t[i - 1];
            dy           = y[i] - y[i - 1];
            tmp_lists(0) = dx;
            for (int j = 1; j < N; j++)
            {
                tmp_lists(j) = tmp_lists(j - 1) * dx;
            }

            A = Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1);
            for (int j = 0; j < N - 1; j++)
            {
                for (int k = j + 1; k < N - 1; k++)
                {
                    A(j, k) = CMatrix(j, k) * tmp_lists(k - j - 1);
                }
                Real tmp = CMatrix(j, N - 1) / tmp_lists(j);
                b(j)     = tmp * dy;
                for (int k = 0; k < N - 1; k++)
                {
                    A(j, k) -= tmp * tmp_lists(k);
                }
            }
            // C_i = A^{-1} C_{i+1} - A^{-1} b
            Eigen::ColPivHouseholderQR<Eigen::Matrix<Real, N - 1, N - 1>> solver;
            solver.compute(A);
            C = solver.solve(C - b);

            tmpx = dx;
            C_iN = dy;
            for (int k = 0; k < N - 1; k++)
            {
                C_iN -= C(k) * tmpx;
                tmpx *= dx;
            }
            C_iN /= tmpx;
            // compute smaller first, and use EIGEN to solve the linear system

            tmp_coeffs[0] = y[i - 1];
            for (int k = 0; k < N - 1; k++)
            {
                tmp_coeffs[k + 1] = C(k);
            }
            tmp_coeffs[N] = C_iN;
            coeffs[i - 1] = tmp_coeffs;
        }

        poly = PPoly<Real>(coeffs, t, 0);
        return;
    }
    else if (method == 1)
    {
        /**
         * @brief all derivatives at t[0] are given in boundary_condition
         */
        int                            t_size = t.size();
        std::vector<std::vector<Real>> coeffs; // coefficients of the polynomial
        Eigen::Matrix<Real, N - 1, 1>  C     = Eigen::Matrix<Real, N - 1, 1>::Zero(N - 1);
        Real                           denom = Real(1.0);
        for (int i = 0; i < N - 1; i++)
        {
            C(i) = boundary_condition[i] / denom;
            denom *= Real(i + 2);
        }
        Real dx   = t[1] - t[0];
        Real dy   = y[1] - y[0];
        Real tmpx = dx;
        Real C_iN = dy;
        for (int i = 0; i < N - 1; i++)
        {
            C_iN -= C(i) * tmpx;
            tmpx *= dx;
        }
        C_iN /= tmpx;
        // coeffs.push_back(C, C_iN);
        std::vector<Real> tmp_coeffs(N, 0.0);
        tmp_coeffs[0] = y[0];
        for (int i = 0; i < N - 1; i++)
        {
            tmp_coeffs[i + 1] = C(i);
        }
        tmp_coeffs.push_back(C_iN);
        coeffs.push_back(tmp_coeffs);
        // first finished, iteration begin.

        Eigen::Matrix<Real, N - 1, N - 1> A = Eigen::Matrix<Real, N - 1, N - 1>::Zero(N - 1, N - 1);
        Eigen::Matrix<Real, N - 1, 1>     b = Eigen::Matrix<Real, N - 1, 1>::Zero(N - 1);
        Eigen::Matrix<Real, N - 1, 1>     tmp_lists = Eigen::Matrix<Real, N - 1, 1>::Ones(N - 1);
        for (int i = 0; i < t_size - 2; i++)
        {
            dx           = t[i + 1] - t[i];
            dy           = y[i + 1] - y[i];
            tmp_lists(0) = dx;
            for (int j = 1; j < N - 1; j++)
            {
                tmp_lists(j) = tmp_lists(j - 1) * dx;
            }

            A = Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1);
            for (int k = 0; k < N - 1; k++)
            {
                for (int j = k + 1; j < N - 1; j++)
                {
                    A(k, j) = CMatrix(k, j) * tmp_lists(j - k - 1); // tmp_lists(j-i) = dx^{j-i}
                }
                // b(i) = CMatrix(i, N - 1) * dy / tmp_lists(i);
                Real tmp = CMatrix(k, N - 1) / tmp_lists(k);
                b(k)     = tmp * dy;
                for (int j = 0; j < N - 1; j++)
                {
                    A(k, j) -= tmp * tmp_lists(j);
                }
            }
            // C    = b + A * C;
            // Eigen::PartialPivLU<Eigen::Matrix<Real, N - 1, N - 1>> solver;
            // solver.compute(Eigen::Matrix<Real, N - 1, N - 1>::Identity(N - 1, N - 1) - A);
            // C    = solver.solve(b);
            C    = b + A * C;
            C_iN = y[i + 2] - y[i + 1];
            dx   = t[i + 2] - t[i + 1];
            tmpx = t[i + 2] - t[i + 1];
            for (int k = 0; k < N - 1; k++)
            {
                C_iN -= C(k) * tmpx;
                tmpx *= dx;
            }
            C_iN /= tmpx;
            tmp_coeffs[N] = C_iN;
            for (int k = 0; k < N - 1; k++)
            {
                tmp_coeffs[k + 1] = C(k);
            }
            tmp_coeffs[0] = y[i + 1];
            coeffs.push_back(tmp_coeffs);
        }

        poly = PPoly<Real>(coeffs, t, 0);
        return;
    }
    return;
}

template <int N, typename Real>
Real
PPInterpolate<N, Real>::operator()(Real x) const
{
    return poly(x);
}

template <int N, typename Real>
PPoly<Real>
PPInterpolate<N, Real>::getPoly() const
{
    return poly;
}
