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
                 const int           n,
                 const int           check)
    : _coeffs(coeffs)
    , _n(n)
{
    _basis.resize(n+1);
    for (int i=0; i<=n; i++)
    {
        _total_basis.push_back(std::vector<double>(i+1, 0.0));
    }

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
    , _basis(other._basis)
    , _total_basis(other._total_basis)
    , _derivative_basis(other._derivative_basis)
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

    // int interval = begin_idx;
    auto it = std::lower_bound(_t.begin() + begin_idx, _t.begin() + end_idx + 1, x);

    // // return the index of the interval
    if (std::abs(*it - _t[begin_idx]) < 1e-10)
    {
        return begin_idx;
    }
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
    
    // ignore the return value of get_basis
    get_basis(x);

    // Use precomputed coefficients directly
    double result = 0.0;
    for (int i=0; i<=_n; i++)
    {
        result += _coeffs[interval - _n + 1 + i] * _basis[i];
    }
    return result;
}


std::vector<double> 
BSpline::get_basis(const double x)
{
    int interval = get_interval(x);
    if (interval == -1)
    {
        // return std::vector<double>(_n + 1, 0.0); // 返回全 0 的向量
        return _basis;
    }


    // std::vector<double> basis(_n + 1, 0.0);
    std::vector<double> left(_n + 1), right(_n + 1);

    // 初始条件
    _basis[0] = 1.0;

    // Calculate basis functions using a single array
    for (int k = 1; k <= _n; ++k)
    {
        for (int j = k; j >= 0; --j) // Reverse order to use previous row values directly
        {
            double left = (j == 0) ? 0 : (x - _t[interval - k + j]) / (_t[interval + j] - _t[interval - k + j]);
            double right = (j == k) ? 0 : (_t[interval + j + 1] - x) / (_t[interval + j + 1] - _t[interval - k + j + 1]);
            
            _basis[j] = left * ((j > 0) ? _basis[j - 1] : 0) + right * ((j < k) ? _basis[j] : 0);
        }
    }

    return _basis;
}



BSpline&
BSpline::operator=(const BSpline &other)
{
    _coeffs = other._coeffs;
    _t = other._t;
    _n = other._n;
    _basis = other._basis;
    _total_basis = other._total_basis;
    _derivative_basis = other._derivative_basis;
    return *this;
}

std::vector<double>
BSpline::get_coeffs() const
{
    return _coeffs;
}

std::vector<double>
BSpline::get_t() const
{
    return _t;
}

int
BSpline::get_n() const
{
    return _n;
}


// std::vector<std::vector<double>>
void
BSpline::get_total_basis(const double x)
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
            double left = (j == 0) ? 0 : (x - _t[interval - k + j]) / (_t[interval + j] - _t[interval - k + j]);
            double right = (j == k) ? 0 : (_t[interval + j + 1] - x) / (_t[interval + j + 1] - _t[interval - k + j + 1]);
            _total_basis[k][j] = left * ((j > 0) ? _total_basis[k - 1][j - 1] : 0) + right * ((j < k) ? _total_basis[k - 1][j] : 0);
            // std::cout << _total_basis[k][j] << std::endl;
        }
    }
}

double
BSpline::cal_derivative_basis(const int interval,
                              const int j,
                              const int degree,
                              const int n)
{
    if (n == 0)
    {
        return _total_basis[degree][j];
    }
    else if (n == 1)
    {
        double left=0.0, right=0.0;
        if (j > 0)
        {
            left = (static_cast<double>(degree) / (_t[interval+j+degree-_n] - _t[interval+j-_n])) * _total_basis[degree-1][j-1];
        }
        if (j < degree)
        {
            right = (static_cast<double>(degree) / (_t[interval+j+degree-_n+1] - _t[interval+j+1-_n]) )* _total_basis[degree-1][j];
        }
        return left - right;
    }
    else
    {
        // double left = static_cast<double>(degree) / (_t[interval+j+degree-_n] - _t[interval+j-_n]) * cal_derivative_basis(interval, j-1, degree-1, n-1);
        // double right = static_cast<double>(degree) / (_t[interval+j+degree-_n+1] - _t[interval+j+1-_n]) * cal_derivative_basis(interval, j, degree-1, n-1);
        double left = 0.0, right=0.0;
        if (j > 0)
        {
            left = (static_cast<double>(degree) / (_t[interval+j+degree-_n] - _t[interval+j-_n])) * cal_derivative_basis(interval, j-1, degree-1, n-1);
        }
        if (j < degree)
        {
            right = (static_cast<double>(degree) / (_t[interval+j+degree-_n+1] - _t[interval+j+1-_n])) * cal_derivative_basis(interval, j, degree-1, n-1);
        }
        return left - right;
    }
}

std::vector<double>
BSpline::basis_derivative(const double x, 
                          const int    n)
{
    int interval = get_interval(x);
    if (interval == -1)
    {
        return std::vector<double>(_n + 1, 0.0);
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


double
BSpline::derivative(const double x,
                    const int    n)
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


    std::vector<double> basis = basis_derivative(x, n);
    double sum = 0.0;
    for (int i = 0; i <= _n; ++i)
    {
        sum += _coeffs[interval - _n + i + 1] * basis[i];
    }
    return sum;
}