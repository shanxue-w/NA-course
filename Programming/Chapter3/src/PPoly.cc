/**
 * @file PPoly.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The implementation of PPoly
 * @version 0.1
 * @date 2024-11-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "PPoly.hpp"
#include <algorithm>
#include <stdexcept>

PPoly::PPoly(){}

PPoly::PPoly(const std::vector<std::vector<double>> &coeffs,
             const std::vector<double> &t,
             const int check)
{
    /**
    * make sure the t is sorted, and the size of t is equal to the size of coeffs.
    * coeffs is corresponding to the interval [t[i], t[i+1]], 
    * so when sorting t, the coeffs should be sorted too.
    * 
    */
    _t = t;
    _coeffs = coeffs;
    if (check)
    {
        // check if the size of t is equal to the size of coeffs
        if (t.size()-1 != coeffs.size())
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

            std::vector<double> new_t(t.size());
            std::vector<std::vector<double>> new_coeffs(t.size());
            for (size_t i = 0; i < t.size(); i++)
            {
                new_t[i] = t[idx[i]];
                new_coeffs[i] = coeffs[idx[i]];
            }
            _t = new_t;
            _coeffs = new_coeffs;
        }
    }
}

PPoly::PPoly(const PPoly &other)
{
    this->_t = other._t;
    this->_coeffs = other._coeffs;
}

int
PPoly::findInterval(double x) const 
{
    /**
     * find the interval that x belongs to.
     * 
     * if x is out of the range of t, return -1.
     * 
     */

    if (x < _t[0] || x > _t[_t.size()-1]) 
        return -1;
    else
    {
        // middle search
        int left = 0;
        int right = _t.size()-1;
        int mid = (left + right) / 2;
        while (right-left > 1)
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


double
PPoly::operator()(double x) const
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
        double xt = x - _t[interval];
        const std::vector<double> &interval_coeffs = _coeffs[interval];
        double result = 0.;
        for (int i=interval_coeffs.size()-1; i>=0; i--)
        {
            result = result * xt + interval_coeffs[i];
        }
        return result;
    }
    return 0.;
}

double _factorial(int n)
{
    if (n == 0)
        return 1.0;
    double result = 1.0;
    for (int i=1; i<=n; i++)
    {
        result *= i;
    }
    return result;
}


double 
PPoly::derivative(double x, int n) const
{
    int idx = findInterval(x);
    int N = _coeffs[idx].size()-1;
    if (n > N)
        return 0.;
    else
    {
        int m = N-n;
        if (idx == -1)
            return 0.;
        else
        {
            std::vector<double> new_t({_t[idx], _t[idx+1]});
            std::vector<std::vector<double>> new_coeffs;
            std::vector<double> new_coeff(m+1);
            double init = _factorial(n); // n!
            for (int i=0; i<m+1; i++)
            {
                new_coeff[i] = _coeffs[idx][i+n] * init;
                init /= (i+1);
                init *= (n+i+1); // init = (i+2)*...*(n+i+1)
            }
            new_coeffs.push_back(new_coeff);
            PPoly new_poly(new_coeffs, new_t);
            return new_poly(x);
        }
    }
    return 0.;
}


double
PPoly::integral(double a, double b) const
{
    /**
     * integrate the polynomial from a to b.
     * 
     * Not implemented yet.
     */
    return (b-a);
}

PPoly &
PPoly::operator=(const PPoly &other)
{
    // construct it 
    if (this != &other)
    {
        this->_t = other._t;
        this->_coeffs = other._coeffs;
    }
    return *this;
}

