/**
 * @file BSpline.tpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief This is the implementation of BSpline
 * @version 0.1
 * @date 2024-11-24
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef BSPLINE_TPP
#define BSPLINE_TPP

#include "BSpline.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

template <typename Real>
BSpline<Real>::BSpline(std::vector<Real> coeffs, std::vector<Real> t, const int n, const int check)
    : _coeffs(coeffs), _n(n)
{
    /**
     * @details For we have only \f$N\f$ knots and \f$N+n-1\f$ coefficients, we
     * need to add \f$2n-2\f$ knots at the beginning and the end of the knot
     *
     * Here simply use equal spaced knots with space \f$1.0\f$
     *
     */
    _basis.resize(n + 1);
    for (int i = 0; i <= n; i++)
    {
        _total_basis.push_back(std::vector<Real>(i + 1, 0.0));
    }

    // Reserve space for _t to prevent multiple reallocations
    _t.reserve(t.size() + 2 * (_n - 1));

    // Prepend and append knots to ensure B-spline boundary conditions
    int len = t.size();
    for (int i = 1 - _n; i < len + (_n - 1); ++i)
    {
        _t.push_back((i < 0) ? t.front() + i : (i >= len) ? t.back() + (i - len + 1) : t[i]);
    }
    len = _t.size();

    // Check if _t.size() and _coeffs.size() align correctly
    if (static_cast<int>(_t.size() - _coeffs.size()) != (_n - 1))
    {
        std::cout << _t.size() << " " << _coeffs.size() << std::endl;
        throw std::invalid_argument("_t.size() - _coeffs.size() != (n-1)");
    }
    // Optional check for sorted knot vector
    if (check == 1 && !std::is_sorted(t.begin(), t.end()))
    {
        throw std::invalid_argument("t must be sorted");
    }
}

template <typename Real>
BSpline<Real>::BSpline(const BSpline &other)
    : _coeffs(other._coeffs), _t(other._t), _basis(other._basis), _total_basis(other._total_basis),
      _derivative_basis(other._derivative_basis), _n(other._n)
{
}

template <typename Real>
int
BSpline<Real>::get_interval(const Real &x) const
{
    /**
     * @details For BSpline, we only care the value in \f$[t_1, t_{N}]\f$,
     * so when x < t_1 or x > t_N, we return -1.
     *
     * Then use `std::lower_bound` to find the first element in _t that is not
     * less than x. Just a binary search.
     */
    int begin_idx = _n - 1;
    int end_idx   = _t.size() - _n;

    // deal with boundary cases
    if (x < _t[begin_idx] || x > _t[end_idx])
    {
        return -1;
    }

    // int interval = begin_idx;
    auto it = std::lower_bound(_t.begin() + begin_idx, _t.begin() + end_idx + 1, x);

    // // return the index of the interval
    if (Real(*it) - Real(_t[begin_idx]) < 1e-10)
    {
        return begin_idx;
    }

    int interval = std::distance(_t.begin(), it) - 1;

    return interval;
}

template <typename Real>
Real
BSpline<Real>::operator()(const Real &x)
{
    /**
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
    int interval = get_interval(x);
    if (interval == -1)
    {
        return 0.0;
    }
    // ignore the return value of get_basis
    get_basis(x);

    // Use precomputed coefficients directly
    Real result = 0.0;
    for (int i = 0; i <= _n; i++)
    {
        result += _coeffs[interval - _n + 1 + i] * _basis[i];
    }
    return result;
}

/**
 * @brief Get the basis functions of the B-spline at the given point, used for `BInterpolate`.
 *
 * @tparam Real the type of the data
 * @param x the point to evaluate the basis functions
 * @return std::vector<Real>
 */
template <typename Real>
std::vector<Real>
BSpline<Real>::get_basis(const Real &x)
{
    int interval = get_interval(x);
    if (interval == -1)
    {
        // return std::vector<Real>(_n + 1, 0.0); // 返回全 0 的向量
        return _basis;
    }

    // std::vector<Real> basis(_n + 1, 0.0);
    std::vector<Real> left(_n + 1), right(_n + 1);

    // 初始条件
    _basis[0] = 1.0;

    // Calculate basis functions using a single array
    for (int k = 1; k <= _n; ++k)
    {
        for (int j = k; j >= 0; --j)
        // Reverse order to use previous row values directly
        {
            Real left = (j == 0) ? Real(0) : (x - _t[interval - k + j]) / (_t[interval + j] - _t[interval - k + j]);
            Real right =
                (j == k) ? Real(0) : (_t[interval + j + 1] - x) / (_t[interval + j + 1] - _t[interval - k + j + 1]);

            _basis[j] = left * ((j > 0) ? _basis[j - 1] : Real(0)) + right * ((j < k) ? _basis[j] : Real(0));
        }
    }

    return _basis;
}

template <typename Real>
BSpline<Real> &
BSpline<Real>::operator=(const BSpline<Real> &other)
{
    _coeffs           = other._coeffs;
    _t                = other._t;
    _n                = other._n;
    _basis            = other._basis;
    _total_basis      = other._total_basis;
    _derivative_basis = other._derivative_basis;
    return *this;
}

template <typename Real>
std::vector<Real>
BSpline<Real>::get_coeffs() const
{
    return _coeffs;
}

/**
 * @brief Get the knot vector of the B-spline
 *
 * @tparam Real the type of the data
 * @return std::vector<Real>
 */
template <typename Real>
std::vector<Real>
BSpline<Real>::get_t() const
{
    return _t;
}

/**
 * @brief Get the degree of the B-spline
 *
 * @tparam Real the type of the data
 * @return int the degree of the B-spline
 */
template <typename Real>
int
BSpline<Real>::get_n() const
{
    return _n;
}

/**
 * @brief Get all available data through the for loop, not just an array
 *
 * @tparam Real the type of the data
 * @param x the point to evaluate
 */
template <typename Real>
void
BSpline<Real>::get_total_basis(const Real &x)
{
    int interval = get_interval(x);
    if (interval == -1)
    {
        return;
    }

    _total_basis[0][0] = 1.0;
    for (int k = 1; k <= _n; ++k)
    {
        for (int j = k; j >= 0; --j)
        {
            Real left = (j == 0) ? Real(0) : (x - _t[interval - k + j]) / (_t[interval + j] - _t[interval - k + j]);
            Real right =
                (j == k) ? Real(0) : (_t[interval + j + 1] - x) / (_t[interval + j + 1] - _t[interval - k + j + 1]);
            _total_basis[k][j] = left * ((j > 0) ? _total_basis[k - 1][j - 1] : Real(0))
                                 + right * ((j < k) ? _total_basis[k - 1][j] : Real(0));
            // std::cout << _total_basis[k][j] << std::endl;
        }
    }
}

/**
 * @brief Calculate the derivative of the basis functions
 *
 * @tparam Real the type of the data
 * @param interval the interval of the point
 * @param j the index of the basis function
 * @param degree the degree of the basis function
 * @param n the order of the derivative
 * @return Real the derivative of the basis function
 */
template <typename Real>
Real
BSpline<Real>::cal_derivative_basis(const int &interval, const int &j, const int &degree, const int &n)
{
    if (n == 0)
    {
        return _total_basis[degree][j];
    }
    else if (n == 1)
    {
        Real left = 0.0, right = 0.0;
        if (j > 0)
        {
            left = (static_cast<Real>(degree) / (_t[interval + j + degree - _n] - _t[interval + j - _n]))
                   * _total_basis[degree - 1][j - 1];
        }
        if (j < degree)
        {
            right = (static_cast<Real>(degree) / (_t[interval + j + degree - _n + 1] - _t[interval + j + 1 - _n]))
                    * _total_basis[degree - 1][j];
        }
        return left - right;
    }
    else
    {
        Real left = 0.0, right = 0.0;
        if (j > 0)
        {
            left = (static_cast<Real>(degree) / (_t[interval + j + degree - _n] - _t[interval + j - _n]))
                   * cal_derivative_basis(interval, j - 1, degree - 1, n - 1);
        }
        if (j < degree)
        {
            right = (static_cast<Real>(degree) / (_t[interval + j + degree - _n + 1] - _t[interval + j + 1 - _n]))
                    * cal_derivative_basis(interval, j, degree - 1, n - 1);
        }
        return left - right;
    }
}

/**
 * @brief Evaluate the derivative of the basis functions at a given point
 *
 * @tparam Real the type of the data
 * @param x the point to evaluate
 * @param n the order of the derivative
 * @return std::vector<Real>
 */
template <typename Real>
std::vector<Real>
BSpline<Real>::basis_derivative(const Real &x, const int &n)
{
    int interval = get_interval(x);
    if (interval == -1)
    {
        return std::vector<Real>(_n + 1, 0.0);
    }
    // std::cout << "interval" << interval << std::endl;

    get_total_basis(x);

    _derivative_basis.resize(_n + 1);

    for (int i = 0; i <= _n; ++i)
    {
        _derivative_basis[i] = cal_derivative_basis(interval, i, _n, n);
        // std::cout << _derivative_basis[i] << std::endl;
    }

    return _derivative_basis;
}

/**
 * @brief Evaluate the derivative of the B-spline at a given point
 *
 * @tparam Real the type of the data
 * @param x the point to evaluate
 * @param n the order of the derivative
 * @return Real the derivative of the B-spline at the given point
 */
template <typename Real>
Real
BSpline<Real>::derivative(const Real &x, const int &n)
{
    if (n > _n)
    {
        return 0.0;
    }
    else if (n < 0)
    {
        throw std::invalid_argument("n must be non-negative");
    }

    int interval = get_interval(x);
    if (interval == -1)
    {
        return 0.0;
    }

    std::vector<Real> basis = basis_derivative(x, n);
    Real              sum   = 0.0;
    for (int i = 0; i <= _n; ++i)
    {
        sum += _coeffs[interval - _n + i + 1] * basis[i];
    }
    return sum;
}

#endif // BSPLINE_TPP