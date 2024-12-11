/**
 * @file PPInterpolate.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implemetnation of PPInterpolate class
 * @version 0.1
 * @date 2024-11-09
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "PPInterpolate.hpp"
#include <Eigen/src/Core/Matrix.h>

/**
 * @brief Construct a new PPInterpolate<1, double>::PPInterpolate object
 *
 * @tparam N The order of the polynomial. Here N = 1.
 * @tparam Real for the type of the coefficients. Here Real = double.
 * @param t the nodes of the interpolation
 * @param y the values of the interpolation
 * @param method methods of interpolation, not used in linear case.
 * @param boundary_condition boundary condition, not used in linear case.
 * @param check whether to check the order of t, default is 0.
 */
template <>
PPInterpolate<1, double>::PPInterpolate(const std::vector<double> &t, // nodes
                                        const std::vector<double> &y, // values
                                        const int                  method,
                                        const std::vector<double> &boundary_condition,
                                        const int                  check) // whether to check the order of t
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    poly = PPoly<double>();
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

/**
 * @brief Interpolate the given data
 *
 * @tparam N for the order of the polynomial. Here N = 1.
 * @tparam Real for the type of the coefficients. Here Real = double.
 * @param t the nodes of the interpolation
 * @param y the values of the interpolation
 * @param method methods of interpolation, not used in linear case.
 * @param boundary_condition boundary condition, not used in linear case.
 */
template <>
void
PPInterpolate<1, double>::interpolate(const std::vector<double> &t, // nodes
                                      const std::vector<double> &y, // values
                                      const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                        // not-a-knot.
                                      const std::vector<double> &boundary_condition) // boundary condition
{
    /**
     * In linear case, we just use the linear function to interpolate the data
     * in each interval.
     *
     * So the form is
     * \f[
     * y = y_i + \frac{y_{i+1} - y_i}{t_{i+1} - t_i} (x - t_i)
     * \f]
     *
     * Store the coefficients and use them to construct the `PPoly` object.
     */
    (void)method;
    (void)boundary_condition;
    std::vector<std::vector<double>> A(t.size() - 1, std::vector<double>(2, 0));

    for (size_t i = 0; i < t.size() - 1; ++i)
    {
        A[i][0] = y[i];
        A[i][1] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }
    poly = PPoly<double>(A, t, 0); // interpolate
}

/**
 * @brief Construct a new PPInterpolate<2, double>::PPInterpolate object
 *
 * @tparam N The order of the polynomial. Here N = 2.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method The method of interpolation, 0 for periodic, 1 for complete.
 * @param boundary_condition The boundary condition.
 * size of the boundary condition should bigger than 0, greater than 1 is ok, but just use the first one.
 * @param check Whether to check the order of t.
 */
template <>
PPInterpolate<2, double>::PPInterpolate(const std::vector<double> &t,      // nodes
                                        const std::vector<double> &y,      // values
                                        const int                  method, // 0 for periodic, 1 for complete.
                                        const std::vector<double> &boundary_condition,
                                        const int                  check) // whether to check the order of t
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    poly = PPoly<double>();
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

/**
 * @brief Interpolate the values of the polynomial at the given nodes.
 *
 * @tparam N The order of the polynomial. Here N = 2.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method The method of interpolation,
 * 0 for periodic, 1 for complete (given one derivative at the beginning point).
 * @param boundary_condition The boundary condition.
 * size of the boundary condition should bigger than 0, greater than 1 is ok, but just use the first one.
 */
template <>
void
PPInterpolate<2, double>::interpolate(const std::vector<double> &t, // nodes
                                      const std::vector<double> &y, // values
                                      const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                        // not-a-knot.
                                      const std::vector<double> &boundary_condition) // boundary condition
{
    /**
     * The Quadratic Spline have the following form:
     * \f[
     * p_i(x) = c_{i,0} + c_{i,1}(x - t_i) + c_{i,2}(x - t_i)^2, x\in [t_i, t_{i+1}].
     * \f]
     *
     * By the values given at the nodes, we have
     * \f[
     * c_{i,0} = p_i(t_i) = y_i,
     * \f]
     *
     * By the continouous of the function, we have
     * \f[
     * y_{i+1} = p_i(t_{i+1}) = y_i + c_{i,1}(t_{i+1} - t_i) + c_{i,2}(t_{i+1} - t_i)^2.
     * \f]
     *
     * To make sure the continuous of the first derivative, we have
     * \f[
     * c_{i+1, 1} = c_{i, 1} + 2c_{i, 2}(t_{i+1} - t_i).
     * \f]
     *
     * Combine the two equations above, we have
     * \f[
     * c_{i,1} + c_{i+1, 1} = 2 f[t_i, t_{i+1}], \quad
     * c_{i, 2} = \frac{ f[t_i, t_{i+1}] - c_{i, 1} } {t_{i+1} - t_i} = \frac{ c_{i+1, 1} - c_{i, 1} } {2 (t_{i+1} - t_i)}.
     * \f]
     *
     * Hence, for \f$c_{i,1}\f$ we already have \f$N-1\f$ equations, and we need one more equation to solve the problem.
     *
     */
    switch (method)
    {
    case 0:
    {
        /**
         * @details **Periodic boundary condition.**
         *
         * The preiodic boundary condition assure that
         * \f[
         * p_1^{\prime}(t_1) = p_{N-1}^{\prime}(t_N) = p_N^{\prime}(t_N) \Rightarrow c_{1,1} = c_{N, 1}.
         * \f]
         *
         * Add the equation to the system, we can solve the coefficients.
         *
         * For the piecewise polynomial only needs \f$ N-1 \f$ polynomials,
         * which explains the loop condition.
         *
         * Finally, we can construct the `PPoly` object.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         *
         */
        int                 t_size = t.size();
        std::vector<double> K_i(t_size - 1, 0);

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        Eigen::SparseMatrix<double> A(t_size, t_size); // t_size items
        Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);

        std::vector<Eigen::Triplet<double>> triplets; // For batch insertion into the sparse matrix

        for (int i = 0; i < t_size - 1; i++)
        {
            b(i) = 2 * K_i[i];

            {
                triplets.push_back(Eigen::Triplet<double>(i, i, 1));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, 1));
            }
        }
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, 0, 1));
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 1, -1));
        b(t_size - 1) = 0;
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        Eigen::VectorXd c = solver.solve(b);

        std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(3, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            coeffs[i][0] = y[i];
            coeffs[i][1] = c(i);
            coeffs[i][2] = (c(i + 1) - c(i)) / (2 * (t[i + 1] - t[i]));
        }
        poly = PPoly<double>(coeffs, t, 0); // interpolate
        return;
    }
    break;

    case 1:
    {
        /**
         * @details **Complete boundary condition.**
         *
         * The complete boundary condition means that all derivatives at t[0] are given in boundary_condition.
         * So \f$c_{1,1} = m_0\f$ is given, directly add it to the system.
         *
         * For the form is very easy, we just use iteration
         * \f[
         * c_{i+1, 1} = 2 f[t_i, t_{i+1}] - c_{i, 1}
         * \f]
         * to solve the coefficients.
         *
         * Finally, we can construct the `PPoly` object.
         */
        int                 t_size = t.size();
        std::vector<double> K_i(t_size - 1, 0);

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        std::vector<double> c(t_size, 0);
        c[0] = boundary_condition[0];

        for (int i = 1; i < t_size; ++i)
        {
            c[i] = 2 * K_i[i - 1] - c[i - 1];
        }

        std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(3, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            coeffs[i][0] = y[i];
            coeffs[i][1] = c[i];
            coeffs[i][2] = (c[i + 1] - c[i]) / (2 * (t[i + 1] - t[i]));
        }
        poly = PPoly<double>(coeffs, t, 0); // interpolate
        return;
    }
    break;

    default:
        /**
         * @throw std::invalid_argument if method is not 0,1
         *
         */
        throw std::invalid_argument("method must be 0,1");
        return;
        break;
    }
    return;
}

/**
 * @brief Construct a new PPInterpolate<3, double>::PPInterpolate object
 *
 * @tparam N The order of the polynomial. Here N = 3.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation.
 * @param y The values of the interpolation.
 * @param method The method used in the interpolation.
 * 0 for periodic, 1 for complete, 2 for natural, others are not supported.
 * @param boundary_condition The boundary condition.
 * size of the boundary condition should bigger than 1, greater than 2 is ok, but just use the first two.
 * @param check Whether to check the order of t.
 */
template <>
PPInterpolate<3, double>::PPInterpolate(const std::vector<double> &t, // nodes
                                        const std::vector<double> &y, // values
                                        const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                          // not-a-knot.
                                        const std::vector<double> &boundary_condition,
                                        const int                  check) // whether to check the order of t
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    poly = PPoly<double>();
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

/**
 * @brief Construct a new PPInterpolate<3, double>::PPInterpolate object
 *
 * @tparam N The order of the polynomial. Here N = 3.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation.
 * @param y The values of the interpolation.
 * @param method The method used in the interpolation.
 * 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
 * @param boundary_condition The boundary condition.
 * size of the boundary condition should bigger than 1, greater than 2 is ok, but just use the first two.
 */
template <>
void
PPInterpolate<3, double>::interpolate(const std::vector<double> &t, // nodes
                                      const std::vector<double> &y, // values
                                      const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                        // not-a-knot.
                                      const std::vector<double> &boundary_condition) // boundary condition
{
    /**
     * The PP-form polynomial have the following form:
     * \f[
     * p_i(x) = c_{i,0} + c_{i,1}(x-x_i) + c_{i,2}(x-x_i)^2 + c_{i,3}(x-x_i)^3, x \in [x_i,x_{i+1}]
     * \f]
     *
     * We let
     * \f[
     * A = \begin{pmatrix}
     * m_1 \\
     * m_2 \\
     * \vdots \\
     * m_n \\
     * \end{pmatrix}
     * \f]
     *
     * We have:
     * \f[
     * \lambda_i m_{i-1} + 2m_i + \mu_i m_{i+1} =
     * 3\mu_i f[x_i, x_{i+1}] + 3\lambda_i f[x_{i-1}, x_i],
     * \f]
     * which gives the first \f$n-2\f$ equations.
     *
     */
    switch (method)
    {
    case 0:
    {
        /**
         *
         * @details **Periodic Cubic Spline**.
         *
         * For periodic boundary condition, we have:
         * \f[
         * m_1 = m_n, \quad M_1 = M_n.
         * \f]
         * The first one is easy to add in the matrix,
         * we only need to represent \f$M_i\f$ by \f$m_i\f$.
         *
         * By simply compute we can have, (\f$M_n\f$ is at the right interval of the last function, so more complex)
         * \f[M_1 = 2\frac{3K_1-2m_1-m_2}{x_2-x_1}, \quad K_1=f[x_1,x_2],\f]
         * \f[M_n = 2\frac{3K_{n-1}-2m_{n-1}-m_n}{x_n-x_{n-1}} + 6\frac{m_{n-1}+m_n-2K_{n-1}}{x_n-x_{n-1}}, \quad K_{n-1} = f[x_{n-1},x_{n}].\f]
         *
         * So we have the last boundary condition as
         * \f[\frac{4}{x_2-x_1}m_1 + \frac{2}{x_2-x_1}m_2 +
         * \frac{2}{x_n-x_{n-1}}m_{n-1} + \frac{4}{x_n-x_{n-1}}m_n =
         * 6\frac{K_1}{x_2-x_1} + 6\frac{K_{n-1}}{x_n-x_{n-1}}.\f]
         *
         * Then we can solve the system of linear equations.
         *
         * Finally compute the coefficients of the cubic spline by
         * \f[
         * c_{i,0} = y_i,\quad
         * c_{i,1} = m_i, \quad
         * c_{i,2} = \frac{3K_i - 2m_{i+1} - m_i}{x_{i+1} - x_i}, \quad
         * c_{i,3} = \frac{m_i - 2K_i + m_{i+1}}{(x_{i+1} - x_i)^2}.
         * \f]
         * Then we can construct the `PPoly` object.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         */

        int                 t_size = t.size();
        std::vector<double> K_i(t_size - 1, 0);

        // Compute K_i = f[x_i, x_{i+1}]

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        Eigen::SparseMatrix<double> A(t_size, t_size);
        Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);

        std::vector<Eigen::Triplet<double>> triplets; // For batch insertion into the sparse matrix
        triplets.reserve(4 * (t_size - 2));

        // Precompute mu_i and lambda_i to avoid redundant calculations

        for (int i = 0; i < t_size - 2; ++i)
        {
            double mu_i     = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            double lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);

            b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);

            // Store the triplets for A matrix in batch
            // Using thread-safe data structure (e.g., critical section or
            // atomic operations) may be required here

            {
                triplets.push_back(Eigen::Triplet<double>(i, i, lambda_i));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, 2.0));
                triplets.push_back(Eigen::Triplet<double>(i, i + 2, mu_i));
            }
        }

        // Boundary condition handling for the last row
        // first derivative to be equal
        triplets.push_back(Eigen::Triplet<double>(t_size - 2, 0, 1.0));
        triplets.push_back(Eigen::Triplet<double>(t_size - 2, t_size - 1, -1.0));
        b(t_size - 2) = 0.0;

        // Boundary condition for the last element in A and b
        // second derivative to be equal
        double dt_0 = t[1] - t[0];
        double dt_n = t[t_size - 1] - t[t_size - 2];
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, 0, 4.0 / dt_0));
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, 1, 2.0 / dt_0));
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 2, 2.0 / dt_n));
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 1, 4.0 / dt_n));

        b(t_size - 1) = (6.0 * K_i[t_size - 2] / dt_n + 6.0 * K_i[0] / dt_0);

        // Fill matrix A in batch
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        Eigen::VectorXd c = solver.solve(b);

        // Prepare the cubic spline coefficients
        std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(4, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            double tmp   = t[i + 1] - t[i];
            coeffs[i][0] = y[i];
            coeffs[i][1] = c(i);
            coeffs[i][2] = (3.0 * K_i[i] - 2.0 * c(i) - c(i + 1)) / tmp;
            coeffs[i][3] = (c(i) - 2.0 * K_i[i] + c(i + 1)) / (tmp * tmp);
        }
        poly = PPoly<double>(coeffs, t, 0); // Interpolate
        return;
    }
    break;

    case 1:
    {
        /**
         *@details **Complete Cubic Spline**
         *
         * Complete Cubic Spline gives the first derivative at the beginning point and
         * the last point.
         * So we have
         * \f[
         * c_{1,1} = m_0, \quad c_{n,1} = m_{n}.
         * \f]
         *
         * Directly add the equations to the system, we can solve the coefficients.
         * Then we can construct the `PPoly` object.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         */
        int                 t_size = t.size();
        std::vector<double> K_i(t_size - 1, 0);

        // compute K_i = f[x_i, x_{i+1}]

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        // init SparseMatrix A and Vector b
        Eigen::SparseMatrix<double>         A(t_size - 2, t_size - 2);
        Eigen::VectorXd                     b = Eigen::VectorXd::Zero(t_size - 2);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(3 * (t_size - 2));

        // 填充矩阵 A 和向量 b

        for (int i = 0; i < t_size - 2; ++i)
        {
            double mu_i     = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            double lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
            if (i == 0)
            {
                b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]) - lambda_i * boundary_condition[0];
            }
            else if (i < t_size - 3)
            {
                b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);
            }
            else
            {
                b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]) - mu_i * boundary_condition[1];
            }

            {
                triplets.emplace_back(i, i, 2.0);
                if (i == 0)
                {
                    triplets.emplace_back(i, i + 1, mu_i);
                }
                else if (i < t_size - 3)
                {
                    triplets.emplace_back(i, i - 1, lambda_i);
                    triplets.emplace_back(i, i + 1, mu_i);
                }
                else
                {
                    triplets.emplace_back(i, i - 1, lambda_i);
                }
            }
        }

        A.setFromTriplets(triplets.begin(), triplets.end());

        // 使用稀疏矩阵求解器
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        Eigen::VectorXd c = solver.solve(b);

        // 构建多项式系数
        std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(4, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            double tmp   = t[i + 1] - t[i];
            coeffs[i][0] = y[i];

            if (i == 0)
            {
                coeffs[i][1] = boundary_condition[0];
                coeffs[i][2] = (3 * K_i[i] - 2 * boundary_condition[0] - c(0)) / tmp;
                coeffs[i][3] = (boundary_condition[0] - 2 * K_i[i] + c(0)) / (tmp * tmp);
            }
            else if (i < t_size - 2)
            {
                coeffs[i][1] = c(i - 1);
                coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - c(i)) / tmp;
                coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + c(i)) / (tmp * tmp);
            }
            else
            {
                coeffs[i][1] = c(i - 1);
                coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - boundary_condition[1]) / tmp;
                coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + boundary_condition[1]) / (tmp * tmp);
            }
        }

        // 创建 PPoly 对象
        poly = PPoly<double>(coeffs, t, 0);
        return;
    }
    break;

    case 2:
    {
        /**
         * @details **Natural Cubic Spline**
         *
         * Natural Cubic Spline means that the second derivative at the beginning point and the last point are zero.
         *
         * Here we use the theorem
         * \f[
         * \mu_i M_{i-1} + 2 M_i + \lambda_i M_{i+1} = 6 f[x_{i-1}, x_{i}, x_{i+1}], \quad i = 2,3,\cdots,n-1.
         * \f]
         *
         * And use
         * \f[
         * M_1 = 0, \quad M_n = 0.
         * \f]
         *
         * We can solve the coefficients.
         *
         * Finally, we compute the coefficients of the pp-form polynomial by
         * \f[
         * c_{i,0} = y_i, \quad c_{i,1} = f[x_i,x_{i+1}] - \frac{1}{6} (M_{i+1}+2M_i)(x_{i+1}-x_i),\quad
         * c_{i,2} = \frac{M_i}{2}, \quad c_{i,3} = \frac{M_{i+1}-M_i}{6(x_{i+1}-x_i)}.
         * \f]
         *
         * For the beginning point and the last point, we make use of
         * \f[
         * M_1 = 0, \quad M_n = 0,
         * \f]
         * to simplify the calculation.
         *
         * Finally we can construct the `PPoly` object.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         *
         */
        int                 t_size = t.size();
        std::vector<double> K_i(t_size - 1, 0);

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }
        std::vector<double> J_i(t_size - 2, 0);

        for (int i = 0; i < t_size - 2; ++i)
        {
            J_i[i] = (K_i[i + 1] - K_i[i]) / (t[i + 2] - t[i]);
        }

        Eigen::VectorXd b = Eigen::VectorXd::Zero(t_size - 2);
        // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(t_size - 2, t_size - 2);
        Eigen::SparseMatrix<double>         A(t_size - 2, t_size - 2);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(3 * (t_size - 2));

        for (int i = 0; i < t_size - 2; ++i)
        {
            b(i) = 6 * J_i[i];

            double mu_i     = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            double lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
            // A(i,i) = 2;

            {
                triplets.emplace_back(i, i, 2.0);
                if (i == 0)
                {
                    // A(i, i+1) = lambda_i;
                    triplets.emplace_back(i, i + 1, lambda_i);
                }
                else if (i < t_size - 3)
                {
                    // A(i, i-1) = mu_i;
                    // A(i, i+1) = lambda_i;
                    triplets.emplace_back(i, i - 1, mu_i);
                    triplets.emplace_back(i, i + 1, lambda_i);
                }
                else
                {
                    // A(i, i-1) = mu_i;
                    triplets.emplace_back(i, i - 1, mu_i);
                }
            }
        }

        A.setFromTriplets(triplets.begin(), triplets.end());
        // solve Ax = b

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        Eigen::VectorXd c = solver.solve(b);

        std::vector<std::vector<double>> coeffs(t.size() - 1, std::vector<double>(4, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            coeffs[i][0] = y[i];
            if (i == 0)
            {

                double tmp   = t[i + 1] - t[i];
                coeffs[i][1] = K_i[i] - 1.0 / 6.0 * c(0) * tmp;
                coeffs[i][2] = 0;
                coeffs[i][3] = c(0) / (6.0 * tmp);
            }
            else if (i < t_size - 2)
            {
                double tmp   = t[i + 1] - t[i];
                coeffs[i][1] = K_i[i] - 1.0 / 6.0 * (c(i) + 2 * c(i - 1)) * tmp;
                coeffs[i][2] = c(i - 1) / 2.0;
                coeffs[i][3] = (c(i) - c(i - 1)) / (6.0 * tmp);
            }
            else
            {
                double tmp   = t[i + 1] - t[i];
                coeffs[i][1] = K_i[i] - 1.0 / 3.0 * c(i - 1) * tmp;
                coeffs[i][2] = c(i - 1) / 2.0;
                coeffs[i][3] = -c(i - 1) / (6.0 * tmp);
            }
        }
        poly = PPoly<double>(coeffs, t, 0); // interpolate
        return;
    }
    break;

    default:
        /**
         * @throw std::invalid_argument if method is not 0,1,2
         *
         */
        throw std::invalid_argument("method must be 0,1,2");
        return;
        break;
    }
    return;
}

template <>
PPInterpolate<1, mpf_class>::PPInterpolate(const std::vector<mpf_class> &t, // nodes
                                           const std::vector<mpf_class> &y, // values
                                           const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                             // not-a-knot.
                                           const std::vector<mpf_class> &boundary_condition,
                                           const int                     check) // whether to check the order of t
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    poly = PPoly<mpf_class>();
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

template <>
void
PPInterpolate<1, mpf_class>::interpolate(const std::vector<mpf_class> &t, // nodes
                                         const std::vector<mpf_class> &y, // values
                                         const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                           // not-a-knot.
                                         const std::vector<mpf_class> &boundary_condition) // boundary condition
{
    // unused variable, how to solve it.
    (void)method;
    (void)boundary_condition;
    std::vector<std::vector<mpf_class>> A(t.size() - 1, std::vector<mpf_class>(2, 0));

    for (size_t i = 0; i < t.size() - 1; ++i)
    {
        A[i][0] = y[i];
        A[i][1] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }
    poly = PPoly<mpf_class>(A, t, 0); // interpolate
}

template <>
PPInterpolate<2, mpf_class>::PPInterpolate(const std::vector<mpf_class> &t, // nodes
                                           const std::vector<mpf_class> &y, // values
                                           const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                             // not-a-knot.
                                           const std::vector<mpf_class> &boundary_condition,
                                           const int                     check) // whether to check the order of t
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    poly = PPoly<mpf_class>();
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

template <>
void
PPInterpolate<2, mpf_class>::interpolate(const std::vector<mpf_class> &t, // nodes
                                         const std::vector<mpf_class> &y, // values
                                         const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                           // not-a-knot.
                                         const std::vector<mpf_class> &boundary_condition) // boundary condition
{
    switch (method)
    {
    case 0:
    {
        int                    t_size = t.size();
        std::vector<mpf_class> K_i(t_size - 1, 0);

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        Eigen::SparseMatrix<mpf_class> A(t_size, t_size);
        // Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size);

        std::vector<Eigen::Triplet<mpf_class>> triplets; // For batch insertion into the sparse matrix

        for (int i = 0; i < t_size - 1; i++)
        {
            b(i) = 2 * K_i[i];

            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, 1));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, 1));
            }
        }
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, 0, 1));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 1, -1));
        b(t_size - 1) = 0;
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd c = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

        std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(3, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            coeffs[i][0] = y[i];
            coeffs[i][1] = c(i);
            coeffs[i][2] = (c(i + 1) - c(i)) / (2 * (t[i + 1] - t[i]));
        }
        poly = PPoly<mpf_class>(coeffs, t, 0); // interpolate
        return;
    }
    break;

    case 1:
    {
        int                    t_size = t.size();
        std::vector<mpf_class> K_i(t_size - 1, 0);

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        std::vector<mpf_class> c(t_size, 0);
        c[0] = boundary_condition[0];

        for (int i = 1; i < t_size; ++i)
        {
            c[i] = 2 * K_i[i - 1] - c[i - 1];
        }

        std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(3, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            coeffs[i][0] = y[i];
            coeffs[i][1] = c[i];
            coeffs[i][2] = (c[i + 1] - c[i]) / (2 * (t[i + 1] - t[i]));
        }
        poly = PPoly<mpf_class>(coeffs, t, 0); // interpolate
        return;
    }
    break;

    default:
        throw std::invalid_argument("method must be 0,1");
        return;
        break;
    }
    return;
}

template <>
PPInterpolate<3, mpf_class>::PPInterpolate(const std::vector<mpf_class> &t, // nodes
                                           const std::vector<mpf_class> &y, // values
                                           const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                             // not-a-knot.
                                           const std::vector<mpf_class> &boundary_condition,
                                           const int                     check) // whether to check the order of t
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    poly = PPoly<mpf_class>();
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

template <>
void
PPInterpolate<3, mpf_class>::interpolate(const std::vector<mpf_class> &t, // nodes
                                         const std::vector<mpf_class> &y, // values
                                         const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                           // not-a-knot.
                                         const std::vector<mpf_class> &boundary_condition) // boundary condition
{
    switch (method)
    {
    case 0:
    {
        int                    t_size = t.size();
        std::vector<mpf_class> K_i(t_size - 1, 0);

        // Compute K_i = f[x_i, x_{i+1}]

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        Eigen::SparseMatrix<mpf_class> A(t_size, t_size);
        // Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size);

        // For batch insertion into the sparse matrix
        std::vector<Eigen::Triplet<mpf_class>> triplets;
        triplets.reserve(4 * (t_size - 2));

        // Precompute mu_i and lambda_i to avoid redundant calculations

        for (int i = 0; i < t_size - 2; ++i)
        {
            mpf_class mu_i     = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            mpf_class lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);

            b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);

            // Store the triplets for A matrix in batch
            // Using thread-safe data structure (e.g., critical section or
            // atomic operations) may be required here

            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, lambda_i));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, 2.0));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, mu_i));
            }
        }

        // Boundary condition handling for the last row
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 2, 0, 1.0));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 2, t_size - 1, -1.0));
        b(t_size - 2) = 0.0;

        // Boundary condition for the last element in A and b
        mpf_class dt_0 = t[1] - t[0];
        mpf_class dt_n = t[t_size - 1] - t[t_size - 2];
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, 0, 4.0 / dt_0));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, 1, 2.0 / dt_0));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 2, 2.0 / dt_n));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 1, 4.0 / dt_n));

        b(t_size - 1) = (6.0 * K_i[t_size - 2] / dt_n + 6.0 * K_i[0] / dt_0);

        // Fill matrix A in batch
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd c = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

        // Prepare the cubic spline coefficients
        std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(4, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            mpf_class tmp = t[i + 1] - t[i];
            coeffs[i][0]  = y[i];
            coeffs[i][1]  = c(i);
            coeffs[i][2]  = (3.0 * K_i[i] - 2.0 * c(i) - c(i + 1)) / tmp;
            coeffs[i][3]  = (c(i) - 2.0 * K_i[i] + c(i + 1)) / (tmp * tmp);
        }
        poly = PPoly<mpf_class>(coeffs, t, 0); // Interpolate
        return;
    }
    break;

    case 1:
    {
        int                    t_size = t.size();
        std::vector<mpf_class> K_i(t_size - 1, 0);

        // compute K_i = f[x_i, x_{i+1}]

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }

        // init SparseMatrix A and Vector b
        Eigen::SparseMatrix<mpf_class> A(t_size - 2, t_size - 2);
        // Eigen::VectorXd                     b = Eigen::VectorXd::Zero(t_size
        // - 2);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size - 2);
        std::vector<Eigen::Triplet<mpf_class>>      triplets;
        triplets.reserve(3 * (t_size - 2));

        // 填充矩阵 A 和向量 b

        for (int i = 0; i < t_size - 2; ++i)
        {
            mpf_class mu_i     = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            mpf_class lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
            b(i)               = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);

            {
                triplets.emplace_back(i, i, 2.0);
                if (i == 0)
                {
                    triplets.emplace_back(i, i + 1, mu_i);
                }
                else if (i < t_size - 3)
                {
                    triplets.emplace_back(i, i - 1, mu_i);
                    triplets.emplace_back(i, i + 1, lambda_i);
                }
                else
                {
                    triplets.emplace_back(i, i - 1, lambda_i);
                }
            }
        }

        A.setFromTriplets(triplets.begin(), triplets.end());

        // 使用稀疏矩阵求解器
        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd c = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

        // 构建多项式系数
        std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(4, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            mpf_class tmp = t[i + 1] - t[i];
            coeffs[i][0]  = y[i];

            if (i == 0)
            {
                coeffs[i][1] = boundary_condition[0];
                coeffs[i][2] = (3 * K_i[i] - 2 * boundary_condition[0] - c(0)) / tmp;
                coeffs[i][3] = (boundary_condition[0] - 2 * K_i[i] + c(0)) / (tmp * tmp);
            }
            else if (i < t_size - 2)
            {
                coeffs[i][1] = c(i - 1);
                coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - c(i)) / tmp;
                coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + c(i)) / (tmp * tmp);
            }
            else
            {
                coeffs[i][1] = c(i - 1);
                coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - boundary_condition[1]) / tmp;
                coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + boundary_condition[1]) / (tmp * tmp);
            }
        }

        // 创建 PPoly 对象
        poly = PPoly<mpf_class>(coeffs, t, 0);
        return;
    }
    break;

    case 2:
    {
        /**
         * Natural Cubic Spline
         */
        int                    t_size = t.size();
        std::vector<mpf_class> K_i(t_size - 1, 0);

        for (int i = 0; i < t_size - 1; ++i)
        {
            K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
        }
        std::vector<mpf_class> J_i(t_size - 2, 0);

        for (int i = 0; i < t_size - 2; ++i)
        {
            J_i[i] = (K_i[i + 1] - K_i[i]) / (t[i + 2] - t[i]);
        }

        // Eigen::VectorXd b = Eigen::VectorXd::Zero(t_size - 2);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size - 2);
        // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(t_size - 2, t_size - 2);
        Eigen::SparseMatrix<mpf_class>         A(t_size - 2, t_size - 2);
        std::vector<Eigen::Triplet<mpf_class>> triplets;
        triplets.reserve(3 * (t_size - 2));

        for (int i = 0; i < t_size - 2; ++i)
        {
            b(i) = 6 * J_i[i];

            mpf_class mu_i     = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            mpf_class lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
            // A(i,i) = 2;

            {
                triplets.emplace_back(i, i, 2.0);
                if (i == 0)
                {
                    // A(i, i+1) = lambda_i;
                    triplets.emplace_back(i, i + 1, lambda_i);
                }
                else if (i < t_size - 3)
                {
                    // A(i, i-1) = mu_i;
                    // A(i, i+1) = lambda_i;
                    triplets.emplace_back(i, i - 1, mu_i);
                    triplets.emplace_back(i, i + 1, lambda_i);
                }
                else
                {
                    // A(i, i-1) = mu_i;
                    triplets.emplace_back(i, i - 1, mu_i);
                }
            }
        }

        A.setFromTriplets(triplets.begin(), triplets.end());
        // solve Ax = b

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd c = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

        std::vector<std::vector<mpf_class>> coeffs(t.size() - 1, std::vector<mpf_class>(4, 0));

        for (int i = 0; i < t_size - 1; ++i)
        {
            coeffs[i][0] = y[i];
            if (i == 0)
            {

                mpf_class tmp = t[i + 1] - t[i];
                coeffs[i][1]  = K_i[i] - 1.0 / 6.0 * c(0) * tmp;
                coeffs[i][2]  = 0;
                coeffs[i][3]  = c(0) / (6.0 * tmp);
            }
            else if (i < t_size - 2)
            {
                mpf_class tmp = t[i + 1] - t[i];
                coeffs[i][1]  = K_i[i] - 1.0 / 6.0 * (c(i) + 2 * c(i - 1)) * tmp;
                coeffs[i][2]  = c(i - 1) / 2.0;
                coeffs[i][3]  = (c(i) - c(i - 1)) / (6.0 * tmp);
            }
            else
            {
                mpf_class tmp = t[i + 1] - t[i];
                coeffs[i][1]  = K_i[i] - 1.0 / 3.0 * c(i - 1) * tmp;
                coeffs[i][2]  = c(i - 1) / 2.0;
                coeffs[i][3]  = -c(i - 1) / (6.0 * tmp);
            }
        }
        poly = PPoly<mpf_class>(coeffs, t, 0); // interpolate
        return;
    }
    break;

    default:
        throw std::invalid_argument("method must be 0,1,2");
        return;
        break;
    }
    return;
}

template class PPInterpolate<1, double>;
template class PPInterpolate<2, double>;
template class PPInterpolate<3, double>;
template class PPInterpolate<4, double>;
template class PPInterpolate<5, double>;
template class PPInterpolate<1, mpf_class>;
template class PPInterpolate<2, mpf_class>;
template class PPInterpolate<3, mpf_class>;
template class PPInterpolate<4, mpf_class>;
template class PPInterpolate<5, mpf_class>;