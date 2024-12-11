/**
 * @file BInterpolate.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief
 * @version 0.1
 * @date 2024-12-08
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "BInterpolate.hpp"
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <gmpxx.h>

/**
 * @brief Construct a new BInterpolate<1, double>::BInterpolate object
 *
 * @tparam N The order of the B-spline. Here N = 1.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method The method of interpolation, not used here in \f$ N = 1 \f$
 * @param boundary_condition The boundary condition, not used here in \f$ N = 1 \f$
 * @param check Whether to check the input data, default is 0
 *
 * @details The constructor of the BInterpolate class, which is used to interpolate
 * the data using B-spline.
 */
template <>
BInterpolate<1, double>::BInterpolate(const std::vector<double> &t,
                                      const std::vector<double> &y,
                                      const int                 &method,
                                      const std::vector<double> &boundary_condition,
                                      const int                  check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<double>();
    /**
     *@brief Make sure the input data is sorted.
     * If not, sort the data.
     */
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
    interpolate(t, y, method, boundary_condition); // begin to interpolate
}

/**
 * @brief Interpolate the data using B-spline.
 *
 * @tparam N The order of the B-splin. Here N = 1.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method The method of interpolation, not used here in \f$ N = 1 \f$
 * @param boundary_condition The boundary condition, not used here in \f$ N = 1 \f$
 *
 * @details In the linear case, we not need to calculte the coefficients, the values are
 * the coefficients. So we just need to construct the `BSpline` object.
 */
template <>
void
BInterpolate<1, double>::interpolate(const std::vector<double> &t,
                                     const std::vector<double> &y,
                                     const int                 &method,
                                     const std::vector<double> &boundary_condition)
{
    (void)method;
    (void)boundary_condition;
    _bspline = BSpline<double>(y, t, 1); // default method is 1
    // std::cout << "Interpolation finished." << std::endl;
}

/**
 * @brief Construct a new BInterpolate<2, double>::BInterpolate object
 *
 * @tparam  N The order of the B-spline. Here N = 2.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method
 * The method of interpolation, 0 for periodic,
 * 1 for complete (given one derivative at the beginning point), 2 for special case.
 * @param boundary_condition The boundary condition, used when method is 1.
 * @param check Whether to check the input data, default is 0.
 */
template <>
BInterpolate<2, double>::BInterpolate(const std::vector<double> &t,
                                      const std::vector<double> &y,
                                      const int                 &method,
                                      const std::vector<double> &boundary_condition,
                                      const int                  check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<double>();
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

/**
 * @brief Construct a new BInterpolate<2, double>::BInterpolate object
 *
 * @tparam N The order of the B-spline. Here N = 2.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t  The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method The method of interpolation, 0 for periodic, 1 for complete, 2 for special case.
 * default is 0.
 * @param boundary_condition The boundary condition, default is {}.
 * size of the boundary condition should bigger than 0, greater than 1 is ok, but just use the first one.
 *
 * @throws std::invalid_argument If the method is not supported.
 */
template <>
void
BInterpolate<2, double>::interpolate(const std::vector<double> &t,
                                     const std::vector<double> &y,
                                     const int                 &method,
                                     const std::vector<double> &boundary_condition)
{
    if (method == 0)
    {
        /**
         *
         * We have \f[ B(x) = \sum_{i=0}^{N} a_i B_{i}^2(x), x\in [t_1, t_N]. \f]
         * So we just need to construct the linear system to solve the coefficients
         * and use them to construct the `BSpline` object.
         *
         * @details **Periodic boundary condition.**
         *
         *\f$ B(x_i) = y_i \f$ is needed, which gives the first t_size equations.
         * Where \f$ B(x_i) = \sum_{j=i-1}^{i+1} c_j B_j(x_i) \f$ (Actually i-1 to i).
         * And the basis functions can be calculated by the `get_basis` function
         * defined in the `BSpline` class.
         *
         * For periodic boundary condition, we need to add one additional equation
         * \f$ B^{\prime} (x_0) = B^{\prime} (x_{t_{size} - 1}) \f$.
         * And we have \f$ B^{\prime} (x) = \sum_{j=i-1}^{i+1} c_j B_j^{\prime}(x), x \in (x_{i}, x_{i+1}] \f$.
         *
         * The primes can be easily calculated by the `basis_derivative` function
         * defined in the `BSpline` class.
         *
         * Finally we have \f$ t_{size} + 1\f$ equations. Then use Eigen to solve the
         * linear system.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         */
        int                                 t_size = t.size();
        Eigen::SparseMatrix<double>         A(t_size + 1, t_size + 1);
        Eigen::VectorXd                     b(t_size + 1);
        std::vector<Eigen::Triplet<double>> triplets;
        std::vector<double>                 tmp_coeffs(t_size + 1, 1.0);
        BSpline<double>                     tmp_spline(tmp_coeffs, t, 2);
        triplets.reserve(3 * (t_size + 1));
        // begin to add B(x_i) = y_i
        for (int i = 0; i < t_size; i++)
        {
            std::vector<double> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
            }
            b(i) = y[i];
        }
        // end of B(x_i) = y_i

        // Begin to add the boundary condition.
        std::vector<double> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);
        std::vector<double> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 1);

        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size, t_size - 1, -diff_basis2[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size, t_size, -diff_basis2[2]));
        b(t_size) = 0;
        // End of adding boundary condition.

        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd     x = solver.solve(b);
        std::vector<double> coeffs(x.data(), x.data() + x.size());
        _bspline = BSpline<double>(coeffs, t, 2);
        return;
    }
    else if (method == 1)
    {
        /**
         *
         * @details **Complete boundary condition.**
         *
         *\f$ B(x_i) = y_i \f$ is needed, which gives the first t_size equations.
         * Where \f$ B(x_i) = \sum_{j=i-1}^{i+1} c_j B_j(x_i) \f$ (Actually i-1 to i).
         * And the basis functions can be calculated by the `get_basis` function
         * defined in the `BSpline` class.
         *
         * For the complete boundary condition, we need to add two additional equations
         * \f[
         * B^{\prime} (x_0) = m_0, \quad B^{\prime} (x_{t_{size} - 1}) = m_{t_{size} - 1}
         * \f]
         * Similarly, we can calculate the derivative of the basis functions by the
         * `basis_derivative` function defined in the `BSpline` class.
         *
         * So we can construct the linear system for the coefficients.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         */
        int                                 t_size = t.size();
        Eigen::SparseMatrix<double>         A(t_size + 1, t_size + 1);
        Eigen::VectorXd                     b(t_size + 1);
        std::vector<Eigen::Triplet<double>> triplets;
        std::vector<double>                 tmp_coeffs(t_size + 1, 1.0);
        BSpline<double>                     tmp_spline(tmp_coeffs, t, 2);
        triplets.reserve(3 * (t_size + 1));
        // begin to add B(x_i) = y_i
        for (int i = 0; i < t_size; i++)
        {
            std::vector<double> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
            }
            b(i) = y[i];
        }
        // end of B(x_i) = y_i

        // Begin to add the boundary condition.
        std::vector<double> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);
        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
        b(t_size) = boundary_condition[0]; // given one at the end point.
        // End of adding boundary condition.

        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd     x = solver.solve(b);
        std::vector<double> coeffs(x.data(), x.data() + x.size());
        _bspline = BSpline<double>(coeffs, t, 2);
        return;
    }
    else if (method == 2)
    {
        /**
         * @details **Special case.**
         *
         * This case is to implement Thm 3.58 to solve question `C`,
         * not too much to say.
         */
        int                 t_size = t.size();
        std::vector<double> new_t(t_size + 1);
        new_t[0] = t[0] - 0.5;
        for (int i = 0; i < t_size; ++i)
        {
            new_t[i + 1] = t[i] + 0.5;
        }
        Eigen::SparseMatrix<double>         A(t_size, t_size);
        Eigen::VectorXd                     b(t_size);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(3 * t_size);
        triplets.push_back(Eigen::Triplet<double>(0, 0, 5.0));
        triplets.push_back(Eigen::Triplet<double>(0, 1, 1.0));
        b(0) = 8.0 * y[0] - 2.0 * boundary_condition[0];
        for (int i = 1; i < t_size - 1; ++i)
        {
            triplets.push_back(Eigen::Triplet<double>(i, i - 1, 1.0));
            triplets.push_back(Eigen::Triplet<double>(i, i, 6.0));
            triplets.push_back(Eigen::Triplet<double>(i, i + 1, 1.0));
            b(i) = 8.0 * y[i];
        }
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 2, 1.0));
        triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 1, 5.0));
        b(t_size - 1) = 8.0 * y[t_size - 1] - 2.0 * boundary_condition[1];
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd     c = solver.solve(b);
        std::vector<double> coeffs(c.data(), c.data() + c.size());
        double              tmp1 = 2.0 * boundary_condition[0] - coeffs[0];
        double              tmp2 = 2.0 * boundary_condition[1] - coeffs[t_size - 1];
        coeffs.insert(coeffs.begin(), tmp1);
        coeffs.push_back(tmp2);
        _bspline = BSpline<double>(coeffs, new_t, 2);
        return;
    }
    else
    {
        /**
         * @details **Other conditions**
         *
         * Just throw an exception.
         *
         */
        throw std::invalid_argument("method must be 0, 1, 2");
    }
    return;
}

/**
 * @brief Construct a new BInterpolate<3, double>::BInterpolate object
 *
 * @tparam N The order of the B-spline. Here N = 3.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method The method of the interpolation.
 * 0 for periodic, 1 for complete, 2 for natural.
 * default is 0.
 * @param boundary_condition The boundary condition. default is {}.
 * size of the boundary condition should bigger than 1, greater than 2 is ok, but just use the first two.
 * @param check Whether to check the input data, default is 0
 */
template <>
BInterpolate<3, double>::BInterpolate(const std::vector<double> &t,
                                      const std::vector<double> &y,
                                      const int                 &method,
                                      const std::vector<double> &boundary_condition,
                                      const int                  check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<double>();
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

/**
 * @brief Construct a new BInterpolate<3, double>::BInterpolate object
 *
 * @tparam N The order of the B-spline. Here N = 3.
 * @tparam Real The type of the coefficients. Here Real = double.
 * @param t The nodes of the interpolation
 * @param y The values of the interpolation
 * @param method The method of the interpolation.
 * 0 for periodic, 1 for complete, 2 for natural.
 * default is 0.
 * @param boundary_condition The boundary condition. default is {}.
 * size of the boundary condition should bigger than 1, greater than 2 is ok, but just use the first two.
 */
template <>
void
BInterpolate<3, double>::interpolate(const std::vector<double> &t,
                                     const std::vector<double> &y,
                                     const int &method, // 0 for periodic, 1 for complete, 2 for natural.
                                     const std::vector<double> &boundary_condition)
{
    if (method == 0)
    {
        /**
         *
         * We have \f[ B(x) = \sum_{i=-1}^{N} a_i B_{i}^3(x), x\in [t_1, t_N]. \f]
         * So we just need to construct the linear system to solve the coefficients
         * and use them to construct the `BSpline` object.
         *
         * @details **Periodic conditions**
         *
         * We have (Actually i-2 to i is enough.) \f[ B(x_i) = \sum_{j=i-2}^{i+1} a_j B_j^3(x_i), \f]
         * So by the value condition \f[ B(x_i) = y_i, \quad i = 1, 2, \cdots, N, \f] we can have the first N equations.
         *
         * More generally, we have (left side is open)
         * \f[ B(x) = \sum_{j=i-2}^{i+1} a_j B_j^3(x), x \in (t_i,t_{i+1}] \f]
         *
         * Due to the implement of the `BSpline` class, in the returned vector of the `get_basis` function,
         * we have that at the beginning point \f$ t_1 \f$, the first three basis are non-zero,
         * and at other points \f$ x \in (t_1, t_N] \f$, the basis of the last three is non-zero.
         * Therefore, in the implementation, we specific the situation of the first point.
         *
         * For the periodic condition, we have
         * \f[ B^{\prime} (x_1) = B^{\prime} (x_N),~ B^{\prime\prime} (x_1) = B^{\prime\prime} (x_N). \f]
         * So we can have the last equations as
         * \f[ \sum_{j=-2}^{1} a_j B_j^{\prime} (x_1) = \sum_{j=N-3}^{N} a_j B_j^{\prime} (x_N), \f]
         * \f[ \sum_{j=-2}^{1} a_j B_j^{\prime\prime} (x_1) = \sum_{j=N-3}^{N} a_j B_j^{\prime\prime} (x_N). \f]
         * Same problem as above, we specific the situation of the first point.
         *
         * Simply use the `basis_derivative` function to get the derivative of the basis functions.
         * we can construct the linear system to solve the coefficients.
         *
         * Finally, we can solve the linear system to get the coefficients.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         *
         */
        int                         t_size = t.size();
        Eigen::SparseMatrix<double> A(t_size + 2, t_size + 2);
        // Eigen::VectorXd b(t_size+2);
        Eigen::Matrix<double, Eigen::Dynamic, 1> b(t_size + 2);
        std::vector<Eigen::Triplet<double>>      triplets;
        triplets.reserve(4 * (t_size + 2));
        std::vector<double> tmp_coeff(t_size + 2, 1.0);
        BSpline<double>     tmp_spline(tmp_coeff, t, 3);
        for (int i = 0; i < t_size; ++i)
        {
            std::vector<double> basis = tmp_spline.get_basis(t[i]);
            // Due to the implement of the `BSpline` class, in the returned vector of the `get_basis` function,
            // we have that at the beginning point \f$ t_1 \f$, the first three basis are non-zero,
            // and at other points \f$ x \in (t_1, t_N] \f$, the basis of the last three is non-zero.
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }
        // first prime condition
        std::vector<double> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);
        std::vector<double> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 1);
        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis1[2]));
        triplets.push_back(Eigen::Triplet<double>(t_size, t_size - 1, -diff_basis2[1])); // negative
        triplets.push_back(Eigen::Triplet<double>(t_size, t_size, -diff_basis2[2]));     // negative
        triplets.push_back(Eigen::Triplet<double>(t_size, t_size + 1, -diff_basis2[3])); // negative
        b(t_size) = 0.0;

        // second prime condition
        diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
        diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, 2, diff_basis1[2]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size - 1, -diff_basis2[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size, -diff_basis2[2]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size + 1, -diff_basis2[3]));
        b(t_size + 1) = 0.0;

        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(b);

        std::vector<double> coeff(x.data(), x.data() + x.size());
        _bspline = BSpline<double>(coeff, t, 3);
        return;
    }
    else if (method == 1)
    {
        /**
         *
         * @details **Complete conditions**
         *
         * The first \f$ N \f$ equations are the same as the periodic case.
         * Just need to add the boundary condition.
         *
         * For the complete boundary condition, we have
         * \f[ B^{\prime} (x_0) = m_0, \quad B^{\prime} (x_{t_{size} - 1}) = m_{t_{size} - 1} \f]
         * So we have
         * \f[
         * \sum_{j=-2}^{1} a_j B_j^{\prime} (x_0) = m_0, \quad
         * \sum_{j=N-3}^{N} a_j B_j^{\prime} (x_{t_{size} - 1}) = m_{t_{size} - 1}
         * \f]
         * Use the `basis_derivative` function to get the derivative of the basis functions,
         * we can construct the linear system to solve the coefficients.
         *
         * Finally, we can solve the linear system to get the coefficients.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         *
         */
        int                         t_size = t.size();
        Eigen::SparseMatrix<double> A(t_size + 2, t_size + 2);
        // Eigen::VectorXd b(t_size+2);
        Eigen::Matrix<double, Eigen::Dynamic, 1> b(t_size + 2);
        std::vector<Eigen::Triplet<double>>      triplets;
        triplets.reserve(4 * (t_size + 2));
        std::vector<double> tmp_coeff(t_size + 2, 1.0);
        BSpline<double>     tmp_spline(tmp_coeff, t, 3);
        for (int i = 0; i < t_size; ++i)
        {
            std::vector<double> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }

        std::vector<double> diff_basis = tmp_spline.basis_derivative(t[0], 1);
        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis[2]));
        b(t_size) = boundary_condition[0];

        diff_basis = tmp_spline.basis_derivative(t[t_size - 1], 1);
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size - 1, diff_basis[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size, diff_basis[2]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size + 1, diff_basis[3]));
        b(t_size + 1) = boundary_condition[1];
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<double>                      coeffs(x.data(), x.data() + x.size());

        // // convert x to std::vector
        // tmp_coeff = std::vector<double>(x.data(), x.data() + x.size());

        // _bspline = BSpline(tmp_coeff, t, 3);
        _bspline = BSpline<double>(coeffs, t, 3);
        return;
    }
    else if (method == 2)
    {
        /**
         *
         * @details **Natural conditions**
         * The first \f$ N \f$ equations are the same as the periodic case.
         * We just need to add two natural conditions.
         * \f[ B^{\prime\prime} (x_0) = 0, \quad B^{\prime\prime} (x_{t_{size} - 1}) = 0 \f]
         * So we have
         * \f[
         * \sum_{j=-2}^{1} a_j B_j^{\prime\prime} (x_0) = 0, \quad
         * \sum_{j=N-3}^{N} a_j B_j^{\prime\prime} (x_{t_{size} - 1}) = 0
         * \f]
         *
         *
         * Finally we can solve the linear system to get the coefficients.
         *
         * @note For the matrix here is a sparse matrix, we use the `Eigen::SparseMatrix` class.
         * And the solver is `Eigen::SparseLU` to solve the sparse linear system.
         */
        int                                      t_size = t.size();
        Eigen::SparseMatrix<double>              A(t_size + 2, t_size + 2);
        Eigen::Matrix<double, Eigen::Dynamic, 1> b(t_size + 2);
        std::vector<Eigen::Triplet<double>>      triplets;
        triplets.reserve(4 * (t_size + 2));
        std::vector<double> tmp_coeff(t_size + 2, 1.0);
        BSpline<double>     tmp_spline(tmp_coeff, t, 3);
        for (int i = 0; i < t_size; ++i)
        {
            std::vector<double> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
                triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }

        std::vector<double> diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
        std::vector<double> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
        // M_0
        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis1[2]));
        b(t_size) = 0.0;
        // M_N
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size - 1, diff_basis2[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size, diff_basis2[2]));
        triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size + 1, diff_basis2[3]));
        b(t_size + 1) = 0.0;

        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(b);

        std::vector<double> coeff(x.data(), x.data() + x.size());
        _bspline = BSpline<double>(coeff, t, 3);
        return;
    }
    else
    {
        throw std::invalid_argument("method must be 0, 1, 2");
    }
}

template <>
BInterpolate<1, mpf_class>::BInterpolate(const std::vector<mpf_class> &t,
                                         const std::vector<mpf_class> &y,
                                         const int                    &method,
                                         const std::vector<mpf_class> &boundary_condition,
                                         const int                     check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<mpf_class>();
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

template <>
void
BInterpolate<1, mpf_class>::interpolate(const std::vector<mpf_class> &t,
                                        const std::vector<mpf_class> &y,
                                        const int                    &method,
                                        const std::vector<mpf_class> &boundary_condition)
{
    (void)method;
    (void)boundary_condition;
    _bspline = BSpline<mpf_class>(y, t, 1); // default method is 1
    // std::cout << "Interpolation finished." << std::endl;
}

template <>
BInterpolate<2, mpf_class>::BInterpolate(const std::vector<mpf_class> &t,
                                         const std::vector<mpf_class> &y,
                                         const int                    &method,
                                         const std::vector<mpf_class> &boundary_condition,
                                         const int                     check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<mpf_class>();
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

template <>
void
BInterpolate<2, mpf_class>::interpolate(const std::vector<mpf_class> &t,
                                        const std::vector<mpf_class> &y,
                                        const int                    &method,
                                        const std::vector<mpf_class> &boundary_condition)
{
    if (method == 0)
    {
        int                                         t_size = t.size();
        Eigen::SparseMatrix<mpf_class>              A(t_size + 1, t_size + 1);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size + 1);
        std::vector<Eigen::Triplet<mpf_class>>      triplets;
        std::vector<mpf_class>                      tmp_coeffs(t_size + 1, 1.0);
        BSpline<mpf_class>                          tmp_spline(tmp_coeffs, t, 2);
        triplets.reserve(3 * (t_size + 1));
        for (int i = 0; i < t_size; i++)
        {
            std::vector<mpf_class> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);
        std::vector<mpf_class> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 1);

        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, t_size - 1, -diff_basis2[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, t_size, -diff_basis2[2]));
        b(t_size) = 0;
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.compute(A);
        // Eigen::VectorXd     x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<mpf_class>                      coeffs(x.data(), x.data() + x.size());
        _bspline = BSpline<mpf_class>(coeffs, t, 2);
        return;

        // periodic, to be done.
    }
    else if (method == 1)
    {
        // given one at the beginning point.
        int                                         t_size = t.size();
        Eigen::SparseMatrix<mpf_class>              A(t_size + 1, t_size + 1);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size + 1);
        std::vector<Eigen::Triplet<mpf_class>>      triplets;
        std::vector<mpf_class>                      tmp_coeffs(t_size + 1, 1.0);
        BSpline<mpf_class>                          tmp_spline(tmp_coeffs, t, 2);
        triplets.reserve(3 * (t_size + 1));
        for (int i = 0; i < t_size; i++)
        {
            std::vector<mpf_class> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);

        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        b(t_size) = boundary_condition[0]; // given one at the end point.
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.compute(A);
        // Eigen::VectorXd     x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<mpf_class>                      coeffs(x.data(), x.data() + x.size());
        _bspline = BSpline<mpf_class>(coeffs, t, 2);
        return;
    }
    else if (method == 2)
    {
        // the example in text book.
        int                    t_size = t.size();
        std::vector<mpf_class> new_t(t_size + 1);
        new_t[0] = t[0] - 0.5;
        for (int i = 0; i < t_size; ++i)
        {
            new_t[i + 1] = t[i] + 0.5;
        }
        Eigen::SparseMatrix<mpf_class>              A(t_size, t_size);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size);
        std::vector<Eigen::Triplet<mpf_class>>      triplets;
        triplets.reserve(3 * t_size);
        triplets.push_back(Eigen::Triplet<mpf_class>(0, 0, 5.0));
        triplets.push_back(Eigen::Triplet<mpf_class>(0, 1, 1.0));
        b(0) = 8.0 * y[0] - 2.0 * boundary_condition[0];
        for (int i = 1; i < t_size - 1; ++i)
        {
            triplets.push_back(Eigen::Triplet<mpf_class>(i, i - 1, 1.0));
            triplets.push_back(Eigen::Triplet<mpf_class>(i, i, 6.0));
            triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, 1.0));
            b(i) = 8.0 * y[i];
        }
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 2, 1.0));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 1, 5.0));
        b(t_size - 1) = 8.0 * y[t_size - 1] - 2.0 * boundary_condition[1];
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.compute(A);
        // Eigen::VectorXd     c = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);
        std::vector<mpf_class>                      coeffs(c.data(), c.data() + c.size());
        mpf_class                                   tmp1 = 2.0 * boundary_condition[0] - coeffs[0];
        mpf_class                                   tmp2 = 2.0 * boundary_condition[1] - coeffs[t_size - 1];
        coeffs.insert(coeffs.begin(), tmp1);
        coeffs.push_back(tmp2);
        _bspline = BSpline<mpf_class>(coeffs, new_t, 2);
        return;
    }
    else
    {
        throw std::invalid_argument("method must be 0, 1, 2");
    }
    return;
}

template <>
BInterpolate<3, mpf_class>::BInterpolate(const std::vector<mpf_class> &t,
                                         const std::vector<mpf_class> &y,
                                         const int                    &method,
                                         const std::vector<mpf_class> &boundary_condition,
                                         const int                     check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<mpf_class>();
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

template <>
void
BInterpolate<3, mpf_class>::interpolate(const std::vector<mpf_class> &t,
                                        const std::vector<mpf_class> &y,
                                        const int &method, // 0 for periodic, 1 for complete, 2 for natural,
                                                           // 3 for not-a-knot.
                                        const std::vector<mpf_class> &boundary_condition)
{
    if (method == 0)
    {
        int                            t_size = t.size();
        Eigen::SparseMatrix<mpf_class> A(t_size + 2, t_size + 2);
        // Eigen::VectorXd b(t_size+2);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size + 2);
        std::vector<Eigen::Triplet<mpf_class>>      triplets;
        triplets.reserve(4 * (t_size + 2));
        std::vector<mpf_class> tmp_coeff(t_size + 2, 1.0);
        BSpline<mpf_class>     tmp_spline(tmp_coeff, t, 3);
        for (int i = 0; i < t_size; ++i)
        {
            std::vector<mpf_class> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }
        std::vector<mpf_class> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);
        std::vector<mpf_class> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 1);
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 2, diff_basis1[2]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, t_size - 1, -diff_basis2[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, t_size, -diff_basis2[2]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, t_size + 1, -diff_basis2[3]));
        b(t_size) = 0.0;

        diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
        diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, 2, diff_basis1[2]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size - 1, -diff_basis2[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size, -diff_basis2[2]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size + 1, -diff_basis2[3]));
        b(t_size + 1) = 0.0;

        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);

        std::vector<mpf_class> coeff(x.data(), x.data() + x.size());
        _bspline = BSpline<mpf_class>(coeff, t, 3);
        return;
    }
    else if (method == 1)
    {
        int                            t_size = t.size();
        Eigen::SparseMatrix<mpf_class> A(t_size + 2, t_size + 2);
        // Eigen::VectorXd b(t_size+2);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size + 2);
        std::vector<Eigen::Triplet<mpf_class>>      triplets;
        triplets.reserve(4 * (t_size + 2));
        std::vector<mpf_class> tmp_coeff(t_size + 2, 1.0);
        BSpline<mpf_class>     tmp_spline(tmp_coeff, t, 3);
        for (int i = 0; i < t_size; ++i)
        {
            std::vector<mpf_class> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis = tmp_spline.basis_derivative(t[0], 1);
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 0, diff_basis[0]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 1, diff_basis[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 2, diff_basis[2]));
        b(t_size) = boundary_condition[0];

        diff_basis = tmp_spline.basis_derivative(t[t_size - 1], 1);
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size - 1, diff_basis[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size, diff_basis[2]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size + 1, diff_basis[3]));
        b(t_size + 1) = boundary_condition[1];
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<mpf_class>                      coeffs(x.data(), x.data() + x.size());

        // // convert x to std::vector
        // tmp_coeff = std::vector<mpf_class>(x.data(), x.data() + x.size());

        // _bspline = BSpline(tmp_coeff, t, 3);
        _bspline = BSpline<mpf_class>(coeffs, t, 3);
        return;
    }
    else if (method == 2)
    {
        int                                         t_size = t.size();
        Eigen::SparseMatrix<mpf_class>              A(t_size + 2, t_size + 2);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size + 2);
        std::vector<Eigen::Triplet<mpf_class>>      triplets;
        triplets.reserve(4 * (t_size + 2));
        std::vector<mpf_class> tmp_coeff(t_size + 2, 1.0);
        BSpline<mpf_class>     tmp_spline(tmp_coeff, t, 3);
        for (int i = 0; i < t_size; ++i)
        {
            std::vector<mpf_class> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[0]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
        std::vector<mpf_class> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 2, diff_basis1[2]));
        b(t_size) = 0.0;
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size - 1, diff_basis2[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size, diff_basis2[2]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size + 1, t_size + 1, diff_basis2[3]));
        b(t_size + 1) = 0.0;

        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);

        std::vector<mpf_class> coeff(x.data(), x.data() + x.size());
        _bspline = BSpline<mpf_class>(coeff, t, 3);
        return;
    }
    else if (method == 3)
    {
        // to be done
    }
    else
    {
        throw std::invalid_argument("method must be 0, 1, 2");
    }
}

template class BInterpolate<1, double>;
template class BInterpolate<2, double>;
template class BInterpolate<3, double>;
template class BInterpolate<1, mpf_class>;
template class BInterpolate<2, mpf_class>;
template class BInterpolate<3, mpf_class>;