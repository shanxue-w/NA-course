/**
 * @file BSpline.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Define the BSpline class.
 * @version 0.1
 * @date 2024-11-11
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef BSPLINE_HPP
#define BSPLINE_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <omp.h>
#include <vector>

#define EIGEN_USE_THREADS
#define EIGEN_USE_LAPACKE

// using Real = double;

/**
 * @brief The definition of the class BSpline, which is a B-spline class.
 *
 * @tparam Real The type of the real number, default is double.
 * Can be double, long double, mpf_class, etc.
 *
 * @details The BSpline is in the form of
 * \f[ B(x) = \sum_{j=i-n+1}^{i+1} c_j B_j^n(x), ~ x\in (x_i, x_{i+1}] \f]
 * In my implementation, the left side is open and the right side is closed.
 * The computation of \f$ B_j^n(x) \f$ is very important.
 * And the main function to do this is `get_basis(const Real &x)`.
 * This function is high performance, which will be explained in the implementation.
 *
 * This class overrides the `operator()` to make it easier to use. Other functions are also provided for convenience.
 */
template <typename Real = double> // template <typename T>, default is double
class BSpline
{
    /**
     * @brief The Bspline is in the form of
     * \f[ B(x) = \sum_{j=i-n+1}^{i+1} c_j B_j^n(x), ~ x\in[x_i, x_{i+1}] \f]
     * So the computation of B_j^n(x) is very important.
     */
private:
    std::vector<Real>              _coeffs;
    std::vector<Real>              _t;
    std::vector<Real>              _basis;
    std::vector<std::vector<Real>> _total_basis;
    std::vector<Real>              _derivative_basis;
    int                            _n;
    std::vector<std::vector<Real>> _difference_t;

public:
    BSpline() = default;

    /**
     * @brief Construct a new BSpline::BSpline object
     *
     * @param coeffs The coefficients of the B-spline
     * @param t \f$[t_1,\cdots,t_{N-1}, t_{N},]\f$.
     * @param n The order of the B-spline
     * @param check Whether to check the input data, default is 0.
     *
     * @throw std::invalid_argument if the size of t is not sorted or the final
     * size of t minus the size of coeffs is not equal to n-1.
     *
     * @details For we have only \f$N\f$ knots and \f$N+n-1\f$ coefficients, we
     * need to add \f$2n-2\f$ knots at the beginning and the end of the knot
     *
     * Here simply use equal spaced knots with space \f$1.0\f$
     *
     */
    BSpline(const std::vector<Real> coeffs, const std::vector<Real> t, const int n, const int check = 0);

    /**
     * @brief Construct a new BSpline object by copying another one.
     *
     * @param other The object to copy, class `BSpline`.
     */
    BSpline(const BSpline &other);

    /**
     * @brief The minimum idx of the interval that contains x
     *
     * @param x The value to be checked
     * @return int The minimum idx of the interval that contains x
     *
     * @details For BSpline, we only care the value in \f$[t_1, t_{N}]\f$,
     * so when \f$x < t_1\f$ or \f$x > t_N\f$, we return -1.
     *
     * Then use `std::lower_bound` to find the first element in _t that is not
     * less than x. Just a binary search.
     */
    int
    get_interval(const Real &x) const;

    /**
     * @brief Overload the assignment operator.
     *
     * @param other The object to copy, class `BSpline`.
     * @return BSpline&
     */
    BSpline &
    operator=(const BSpline &other);

    /**
     * @brief Overload the operator `()` to evaluate the B-spline at x.
     * @param x The value to be evaluated
     * @return Real The value of the B-spline at x
     */
    Real
    operator()(const Real &x);

    /**
     * @brief Get the coeffs object
     *
     * @return std::vector<Real>
     *
     * @details For every point \f$x \in(x_i, x_{i+1}]\f$, we know there are only \f$n+1\f$ basis functions are non-zero.
     * More specifically, they are
     * \f[
     * B_{i-n+1}^n, B_{i-n+2}^n, \cdots, B_{i+1}^n.
     * \f]
     * Therefore, we only need to calculate the \f$n+1\f$ basis functions and multiply them with the corresponding coefficients,
     * then sum them up we can get the value of the B-spline at \f$x\f$.
     *
     * Calculate the basis functions using a single array and use the formula
     * \f[
     * B_i^n(x) = \frac{x - t_{i-1}}{t_{i+n-1} - t_{i-1}} B_{i}^{n-1}(x) + \frac{t_{i+n} - x}{t_{i+n} - t_{i}} B_{i+1}^{n-1}(x).
     * \f]
     *
     */
    std::vector<Real>
    get_coeffs() const;

    /**
     * @brief Get the t object
     *
     * @return std::vector<Real>
     */
    std::vector<Real>
    get_t() const;

    /**
     * @brief Get the n object, which is the order of the B-spline
     * @return int
     */
    int
    get_n() const;

    /**
     * @brief Get the basis developed by a signal point x, which will disscussed in `operator()` function.
     *
     * @param x The signal point
     * @return std::vector<Real>
     */
    std::vector<Real>
    get_basis(const Real &x);

    /**
     * @brief Get the total basis object.
     *
     * @param x The signal point
     */
    void
    get_total_basis(const Real &x);

    /**
     * @brief Get the derivative basis object.
     *
     * @param interval The interval that contains x
     * @param j The index of the basis function
     * @param degree The order of the B-spline
     * @param n The index of the derivative
     * @return Real The value of the derivative basis function
     */
    Real
    cal_derivative_basis(const int &interval, const int &j, const int &degree, const int &n);

    /**
     * @brief Get the derivative basis object.
     *
     * @param x The signal point
     * @param n The index of the derivative
     * @return std::vector<Real>
     */
    std::vector<Real>
    basis_derivative(const Real &x, const int &n);

    /**
     * @brief Evaluate the n-th derivative of the B-spline at x.
     *
     * @param x The value to be evaluated
     * @param n The order of the derivative
     * @return Real The value of the n-th derivative of the B-spline at x
     */
    Real
    derivative(const Real &x, const int &n = 1);
};

#pragma once
#include "BSpline.tpp"

#endif // BSPLINE_HPP