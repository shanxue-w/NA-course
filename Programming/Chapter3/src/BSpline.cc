/**
 * @file BSpline.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implementation of BSpline class
 * @version 0.1
 * @date 2024-11-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "BSpline.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>


/**
 * @brief Construct a new BSpline::BSpline object
 * 
 * @param coeffs 
 * @param t \f$[t_0,\cdots,t_{n-2}, t_{n-1},\cdots,t_{len-1}]\f$, the first n-1 elements are added, we actually interpolate at f(t_{n-1}), \cdots, f(t_{len-1})
 * @param n 
 * @param check 
 */
BSpline::BSpline(std::vector<double> coeffs,
                 std::vector<double> t,
                 const int n,
                 const int check)
    : _coeffs(std::move(coeffs))
    , _n(n)
{
    // Reserve space for _t to prevent multiple reallocations
    _t.reserve(t.size() + 2 * (_n - 1));

    // Prepend and append knots to ensure B-spline boundary conditions
    int len = t.size();
    for (int i = 1 - _n; i < len + (_n - 1); ++i)
    {
        _t.push_back((i < 0) ? t.front() + i : (i >= len) ? t.back() + (i - len + 1) : t[i]);
    }

    // Check if _t.size() and _coeffs.size() align correctly
    if (static_cast<int>(_t.size() - _coeffs.size()) != (_n - 1))
    {
        throw std::invalid_argument("_t.size() - _coeffs.size() != (n-1)");
    }

    // Optional check for sorted knot vector
    if (check == 1 && !std::is_sorted(t.begin(), t.end()))
    {
        throw std::invalid_argument("t must be sorted");
    }
}


BSpline::BSpline(const BSpline &other)
    : _coeffs(other._coeffs)
    , _t(other._t)
    , _n(other._n)
{}

/**
 * @brief The minimum idx of the interval that contains x
 * 
 * @param x 
 * @return int 
 */
int 
BSpline::get_interval(const double x) const
{
    int begin_idx = _n - 1;
    int end_idx = _t.size() - _n;

    // deal with boundary cases
    if (x < _t[begin_idx] || x > _t[end_idx]) 
    {
        return -1;
    }

    // use binary search to find the interval, std::lower_bound is used here
    auto it = std::lower_bound(_t.begin() + begin_idx, _t.begin() + end_idx + 1, x);

    // return the index of the interval
    int interval = std::distance(_t.begin(), it) - 1;
    return interval;
}

double 
BSpline::operator()(const double x)
{
    int interval = get_interval(x);
    if (interval == -1)
    {
        return 0.0;
    }

    // Precomputed storage to avoid repeated allocations
    static thread_local Eigen::VectorXd basis(_n + 1);
    basis.setZero();
    basis[0] = 1.0;

    // Calculate basis functions using a single array
    for (int k = 1; k <= _n; ++k)
    {
        for (int j = k; j >= 0; --j) // Reverse order to use previous row values directly
        {
            double left = (j == 0) ? 0 : (x - _t[interval - k + j]) / (_t[interval + j] - _t[interval - k + j]);
            double right = (j == k) ? 0 : (_t[interval + j + 1] - x) / (_t[interval + j + 1] - _t[interval - k + j + 1]);
            
            basis[j] = left * ((j > 0) ? basis[j - 1] : 0) + right * basis[j];
        }
    }

    // Use precomputed coefficients directly
    Eigen::Map<Eigen::VectorXd> coeffs(&_coeffs[interval - _n + 1], _n + 1);
    double result = coeffs.dot(basis);

    return result;
}

