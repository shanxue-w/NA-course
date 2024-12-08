#include "BInterpolate.hpp"
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <gmpxx.h>

template <>
BInterpolate<1, double>::BInterpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int                 &method,
    const std::vector<double> &boundary_condition,
    const int                  check)
    : _t(t), _y(y), _method(method), _boundary_condition(boundary_condition)
{
    _bspline = BSpline<double>();
    /**
     *@brief Make sure the input data is sorted.
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

template <>
void
BInterpolate<1, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int                 &method,
    const std::vector<double> &boundary_condition)
{
    (void)method;
    (void)boundary_condition;
    _bspline = BSpline<double>(y, t, 1); // default method is 1
    // std::cout << "Interpolation finished." << std::endl;
}

template <>
BInterpolate<2, double>::BInterpolate(
    const std::vector<double> &t,
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

template <>
void
BInterpolate<2, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int                 &method,
    const std::vector<double> &boundary_condition)
{
    if (method == 0)
    {
        /**
         *@brief Periodic boundary condition.
         *
         *@details \f$ B(x_i) = y_i \f$ is needed, which gives the first
         * t_size equations.
         * For periodic boundary condition, we need to add one additional equation
         * \f$ B^{\prime} (x_0) = B^{\prime} (x_{t_size - 1}) \f$.
         *
         * Finally we have t_size + 1 equations. Then use Eigen to solve the
         * linear system.
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
        std::vector<double> diff_basis2 =
            tmp_spline.basis_derivative(t[t_size - 1], 1);

        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size, t_size - 1, -diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size, t_size, -diff_basis2[2]));
        b(t_size) = 0;
        // End of adding boundary condition.

        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd     x = solver.solve(b);
        std::vector<double> coeffs(x.data(), x.data() + x.size());
        //< Construct the BSpline object. >/
        _bspline = BSpline<double>(coeffs, t, 2);
        return;
    }
    else if (method == 1)
    {
        /**
         *@brief The complete boundary condition, give one derivative at the
         * beginning point.
         *
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
         *@brief The theorem in the text book, special case, for problem C and D.
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
        double tmp2 = 2.0 * boundary_condition[1] - coeffs[t_size - 1];
        coeffs.insert(coeffs.begin(), tmp1);
        coeffs.push_back(tmp2);
        _bspline = BSpline<double>(coeffs, new_t, 2);
        return;
    }
    else
    {
        throw std::invalid_argument("method must be 0, 1, 2");
    }
    return;
}

template <>
BInterpolate<3, double>::BInterpolate(
    const std::vector<double> &t,
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

template <>
void
BInterpolate<3, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int &method, // 0 for periodic, 1 for complete, 2 for natural,
                       // 3 for not-a-knot.
    const std::vector<double> &boundary_condition)
{
    if (method == 0)
    {
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
        std::vector<double> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);
        std::vector<double> diff_basis2 =
            tmp_spline.basis_derivative(t[t_size - 1], 1);
        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis1[2]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size, t_size - 1, -diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size, t_size, -diff_basis2[2]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size, t_size + 1, -diff_basis2[3]));
        b(t_size) = 0.0;

        diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
        diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, 0, diff_basis1[0]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, 1, diff_basis1[1]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, 2, diff_basis1[2]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size - 1, -diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size, -diff_basis2[2]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size + 1, -diff_basis2[3]));
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
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size - 1, diff_basis[1]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size, diff_basis[2]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size + 1, diff_basis[3]));
        b(t_size + 1) = boundary_condition[1];
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<double> coeffs(x.data(), x.data() + x.size());

        // // convert x to std::vector
        // tmp_coeff = std::vector<double>(x.data(), x.data() + x.size());

        // _bspline = BSpline(tmp_coeff, t, 3);
        _bspline = BSpline<double>(coeffs, t, 3);
        return;
    }
    else if (method == 2)
    {
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
        std::vector<double> diff_basis2 =
            tmp_spline.basis_derivative(t[t_size - 1], 2);
        triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
        triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis1[2]));
        b(t_size) = 0.0;
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size - 1, diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size, diff_basis2[2]));
        triplets.push_back(
            Eigen::Triplet<double>(t_size + 1, t_size + 1, diff_basis2[3]));
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
    else if (method == 3)
    {
        // to be done
    }
    else
    {
        throw std::invalid_argument("method must be 0, 1, 2");
    }
}

template <>
BInterpolate<1, mpf_class>::BInterpolate(
    const std::vector<mpf_class> &t,
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

template <>
void
BInterpolate<1, mpf_class>::interpolate(
    const std::vector<mpf_class> &t,
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
BInterpolate<2, mpf_class>::BInterpolate(
    const std::vector<mpf_class> &t,
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

template <>
void
BInterpolate<2, mpf_class>::interpolate(
    const std::vector<mpf_class> &t,
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
        BSpline<mpf_class> tmp_spline(tmp_coeffs, t, 2);
        triplets.reserve(3 * (t_size + 1));
        for (int i = 0; i < t_size; i++)
        {
            std::vector<mpf_class> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[0]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis1 =
            tmp_spline.basis_derivative(t[0], 1);
        std::vector<mpf_class> diff_basis2 =
            tmp_spline.basis_derivative(t[t_size - 1], 1);

        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, t_size - 1, -diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, t_size, -diff_basis2[2]));
        b(t_size) = 0;
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.compute(A);
        // Eigen::VectorXd     x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<mpf_class> coeffs(x.data(), x.data() + x.size());
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
        BSpline<mpf_class> tmp_spline(tmp_coeffs, t, 2);
        triplets.reserve(3 * (t_size + 1));
        for (int i = 0; i < t_size; i++)
        {
            std::vector<mpf_class> basis = tmp_spline.get_basis(t[i]);
            if (i == 0)
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[0]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis1 =
            tmp_spline.basis_derivative(t[0], 1);

        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        b(t_size) = boundary_condition[0]; // given one at the end point.
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.compute(A);
        // Eigen::VectorXd     x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<mpf_class> coeffs(x.data(), x.data() + x.size());
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
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size - 1, t_size - 2, 1.0));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size - 1, t_size - 1, 5.0));
        b(t_size - 1) = 8.0 * y[t_size - 1] - 2.0 * boundary_condition[1];
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.compute(A);
        // Eigen::VectorXd     c = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);
        std::vector<mpf_class> coeffs(c.data(), c.data() + c.size());
        mpf_class              tmp1 = 2.0 * boundary_condition[0] - coeffs[0];
        mpf_class tmp2 = 2.0 * boundary_condition[1] - coeffs[t_size - 1];
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
BInterpolate<3, mpf_class>::BInterpolate(
    const std::vector<mpf_class> &t,
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

template <>
void
BInterpolate<3, mpf_class>::interpolate(
    const std::vector<mpf_class> &t,
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
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }
        std::vector<mpf_class> diff_basis1 =
            tmp_spline.basis_derivative(t[0], 1);
        std::vector<mpf_class> diff_basis2 =
            tmp_spline.basis_derivative(t[t_size - 1], 1);
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 2, diff_basis1[2]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, t_size - 1, -diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, t_size, -diff_basis2[2]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, t_size + 1, -diff_basis2[3]));
        b(t_size) = 0.0;

        diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
        diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, 0, diff_basis1[0]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, 1, diff_basis1[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, 2, diff_basis1[2]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size - 1, -diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size, -diff_basis2[2]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size + 1, -diff_basis2[3]));
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
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis =
            tmp_spline.basis_derivative(t[0], 1);
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 0, diff_basis[0]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 1, diff_basis[1]));
        triplets.push_back(Eigen::Triplet<mpf_class>(t_size, 2, diff_basis[2]));
        b(t_size) = boundary_condition[0];

        diff_basis = tmp_spline.basis_derivative(t[t_size - 1], 1);
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size - 1, diff_basis[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size, diff_basis[2]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size + 1, diff_basis[3]));
        b(t_size + 1) = boundary_condition[1];
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        // Eigen::VectorXd x = solver.solve(b);
        Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> x = solver.solve(b);
        std::vector<mpf_class> coeffs(x.data(), x.data() + x.size());

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
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 2, basis[2]));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<mpf_class>(i, i, basis[1]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 1, basis[2]));
                triplets.push_back(
                    Eigen::Triplet<mpf_class>(i, i + 2, basis[3]));
            }
            b(i) = y[i];
        }

        std::vector<mpf_class> diff_basis1 =
            tmp_spline.basis_derivative(t[0], 2);
        std::vector<mpf_class> diff_basis2 =
            tmp_spline.basis_derivative(t[t_size - 1], 2);
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 0, diff_basis1[0]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 1, diff_basis1[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size, 2, diff_basis1[2]));
        b(t_size) = 0.0;
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size - 1, diff_basis2[1]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size, diff_basis2[2]));
        triplets.push_back(
            Eigen::Triplet<mpf_class>(t_size + 1, t_size + 1, diff_basis2[3]));
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