/**
 * @file PPoly.tpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implementation of the class PPoly
 * @version 0.1
 * @date 2024-11-24
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef PPOLY_TPP
#define PPOLY_TPP

#include "PPoly.hpp"
#include <algorithm>
#include <stdexcept>

/**
 * @brief Construct a new PPoly<Real>::PPoly object
 *
 * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
 *
 * @details The default constructor of the class PPoly.
 * The coefficients and knots are initialized as empty.
 */
template <typename Real>
PPoly<Real>::PPoly()
{
    _t.clear();
    _coeffs.clear();
}

/**
 * @brief Construct a new PPoly<Real>::PPoly object
 *
 * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
 * @param coeffs The coefficients of pp-form
 * @param t The knots of pp-form
 * @param check Whether to check the input data, default is 0.
 *
 * @details The pp-form polynomial has the form
 * \f[
 * p_i(x) = \sum_{j=0}^{n} a_{i,j} (x - t_i)^j, x \in [t_i, t_{i+1}].
 * \f]
 * \f$a_{i,j}\f$ are the `coeffs`, and \f$t_i\f$ are the `t`.
 */
template <typename Real>
PPoly<Real>::PPoly(const std::vector<std::vector<Real>> &coeffs, const std::vector<Real> &t, const int check)
{
    /**
     * make sure the t is sorted, and the size of t is equal to the size of
     * coeffs. coeffs is corresponding to the interval [t[i], t[i+1]], so when
     * sorting t, the coeffs should be sorted too.
     */
    _t      = t;
    _coeffs = coeffs;
    if (check)
    {
        // check if the size of t is equal to the size of coeffs
        if (t.size() - 1 != coeffs.size())
        {
            throw std::invalid_argument("The size of t is not equal to the size of coeffs");
        }
        // check if t is sorted
        if (!std::is_sorted(t.begin(), t.end()))
        {
            std::vector<size_t> idx(t.size());
            for (size_t i = 0; i < t.size(); i++)
            {
                idx[i] = i;
            }

            std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) { return t[i] < t[j]; });

            for (size_t i = 0; i < t.size(); i++)
            {
                _t[i]      = t[idx[i]];
                _coeffs[i] = coeffs[idx[i]];
            }
        }
    }
}

/**
 * @brief Construct a new PPoly<Real>::PPoly object
 *
 * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
 * @param other The object to copy, class `PPoly`.
 */
template <typename Real>
PPoly<Real>::PPoly(const PPoly<Real> &other)
{
    this->_t      = other._t;
    this->_coeffs = other._coeffs;
}

/**
 * @brief Find the interval that x belongs to.
 *
 * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
 * @param x The value to find the interval.
 * @return int
 */
template <typename Real>
int
PPoly<Real>::findInterval(Real x) const
{
    /**
     *
     * @details When we input a value \f$ x \f$, we need to find the interval that
     * \f$ x \in[t_i, t_{i+1}) \f$, and return the index \f$ i \f$.
     * For the end point \f$ t_n \f$, we return \f$ n-1 \f$.
     *
     * If x is out of the range of the knots, return -1, means
     * we don't have the interval that x belongs to.
     *
     * @note For the knots already sorted, we can use binary search to find the
     * interval.
     *
     */
    if (x < _t[0] || x > _t[_t.size() - 1])
        return -1;
    else
    {
        // middle search
        int left  = 0;
        int right = _t.size() - 1;
        int mid   = (left + right) / 2;
        while (right - left > 1)
        {
            if (x <= _t[mid])
                right = mid;
            else if (x >= _t[mid])
                left = mid;
            else
                return mid;
            mid = (left + right) / 2;
        }
        return mid;
    }
    return -1;
}

/**
 * @brief Overload the operator() to calculate the value of the polynomial at x.
 *
 * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
 * @param x The value to calculate the polynomial.
 * @return Real The value of the polynomial at x.
 */
template <typename Real>
Real
PPoly<Real>::operator()(Real x) const
{
    /**
     *
     * @details The function return the value of the polynomial at a given point.
     *
     * First find the interval that x belongs to.
     * Then use the Horner's method to calculate the value of the polynomial.
     *
     * \f[
     * p_i(x) = \sum_{j=0}^{n} a_{i,j} (x - t_i)^j
     * = a_{i,0} + (x-t_i)(a_{i,1} + (x-t_i)(a_{i,2} + \cdots ))
     * \f]
     */
    int interval = findInterval(x);
    if (interval == -1)
        return 0.;
    else
    {
        Real xt = x - _t[interval];
        // const std::vector<Real> &interval_coeffs = _coeffs[interval];
        // Eigen::VectorXd interval_coeffs = Eigen::Map<const Eigen::VectorXd>(
        //     _coeffs[interval].data(), _coeffs[interval].size());
        Real result = 0.;
        for (int i = _coeffs[interval].size() - 1; i >= 0; i--)
        {
            result = result * xt + _coeffs[interval][i];
        }
        // for (int i = interval_coeffs.size() - 1; i >= 0; i--) {
        //   result = result * xt + interval_coeffs[i];
        // }
        return result;
    }
    return 0.;
}

/**
 * @brief Calculate the n-th derivative of the polynomial at x.
 *
 * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
 * @param x The value to calculate the derivative.
 * @param n The order of the derivative.
 * @return Real The value of the n-th derivative of the polynomial at x.
 *
 * @details If the derivative order n is greater than the degree of the polynomial,
 * return 0. If n is negative, return 0. If x is out of bounds, return 0.
 *
 * Then just derive the polynomial and use Horner's method to evaluate the n-th
 * derivative at x.
 */
template <typename Real>
Real
PPoly<Real>::derivative(Real x, int n) const
{
    // Find the interval that contains x
    int idx = findInterval(x);
    if (idx == -1)
        return 0.0; // If x is out of bounds, return 0

    int N = _coeffs[idx].size() - 1; // Degree of the polynomial in this interval
    if (n > N || n < 0)
        return 0.0; // If the derivative order n is greater than the degree,
                    // return 0

    // The number of coefficients after taking the n-th derivative
    int m = N - n;

    Real xt     = x - _t[idx]; // Calculate (x - t[idx]) for Horner's method
    Real result = 0.0;

    // Initialize factorial multiplier for the n-th derivative
    Real init = 1.0;
    for (Real i = m + 1; i <= N; ++i)
    {
        init *= i; // Compute (N-n+1) * (N-n+2) * ... * N
    }

    // Eigen::VectorXd tmp = Eigen::Map<const Eigen::VectorXd>(&_coeffs[idx][n],
    // m
    // + 1);
    Eigen::Matrix<Real, Eigen::Dynamic, 1> tmp =
        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>>(&_coeffs[idx][n], m + 1);

    // Use Horner's method to evaluate the n-th derivative at x
    for (int i = m; i >= 0; --i)
    {
        // Evaluate polynomial derivative using Horner's method
        result = result * xt + tmp(i) * init;

        // Update the factorial part for the next coefficient
        init /= (n + i);
        init *= i; // Update the factorial part for the next coefficient
    }
    return result;
}

// /**
//  * @brief
//  *
//  * @tparam Real
//  * @param a
//  * @param b
//  * @return Real
//  */
// template <typename Real>
// Real
// PPoly<Real>::integral(Real a, Real b) const
// {
//     /**
//      * integrate the polynomial from a to b.
//      *
//      * Not implemented yet.
//      */
//     return (b - a);
// }

/**
 * @brief Overload the operator= to copy the object.
 *
 * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
 * @param other The object to copy.
 * @return PPoly<Real>&
 */
template <typename Real>
PPoly<Real> &
PPoly<Real>::operator=(const PPoly<Real> &other)
{
    // construct it
    if (this != &other)
    {
        this->_t      = other._t;
        this->_coeffs = other._coeffs;
    }
    return *this;
}

#endif