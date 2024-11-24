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

#include "PPoly.hpp"
#include <algorithm>
#include <stdexcept>

template <typename Real>
PPoly<Real>::PPoly()
{
}

template <typename Real>
PPoly<Real>::PPoly(const std::vector<std::vector<Real>> &coeffs,
                   const std::vector<Real>              &t,
                   const int                             check)
{
    /**
     * make sure the t is sorted, and the size of t is equal to the size of
     * coeffs. coeffs is corresponding to the interval [t[i], t[i+1]], so when
     * sorting t, the coeffs should be sorted too.
     *
     */
    _t      = t;
    _coeffs = coeffs;
    if (check)
    {
        // check if the size of t is equal to the size of coeffs
        if (t.size() - 1 != coeffs.size())
        {
            throw std::invalid_argument(
                "The size of t is not equal to the size of coeffs");
        }
        // check if t is sorted
        if (!std::is_sorted(t.begin(), t.end()))
        {
            std::vector<size_t> idx(t.size());
            for (size_t i = 0; i < t.size(); i++)
            {
                idx[i] = i;
            }

            std::sort(idx.begin(),
                      idx.end(),
                      [&](size_t i, size_t j) { return t[i] < t[j]; });

            for (size_t i = 0; i < t.size(); i++)
            {
                _t[i]      = t[idx[i]];
                _coeffs[i] = coeffs[idx[i]];
            }
        }
    }
}

template <typename Real>
PPoly<Real>::PPoly(const PPoly<Real> &other)
{
    this->_t      = other._t;
    this->_coeffs = other._coeffs;
}

template <typename Real>
int
PPoly<Real>::findInterval(Real x) const
{
    /**
     * find the interval that x belongs to.
     *
     * if x is out of the range of t, return -1.
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

template <typename Real>
Real
PPoly<Real>::operator()(Real x) const
{
    /**
     * First find the interval that x belongs to.
     * Then use the Horner's method to calculate the value of the polynomial.
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

template <typename Real>
Real
PPoly<Real>::derivative(Real x, int n) const
{
    // Find the interval that contains x
    int idx = findInterval(x);
    if (idx == -1)
        return 0.0; // If x is out of bounds, return 0

    int N =
        _coeffs[idx].size() - 1; // Degree of the polynomial in this interval
    if (n > N)
        return 0.0; // If the derivative order n is greater than the degree,
                    // return
                    // 0

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
        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>>(
            &_coeffs[idx][n], m + 1);

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

template <typename Real>
Real
PPoly<Real>::integral(Real a, Real b) const
{
    /**
     * integrate the polynomial from a to b.
     *
     * Not implemented yet.
     */
    return (b - a);
}

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