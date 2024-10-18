#include "Hermite.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include "Polynomial.hpp"

template <class Poly>
Hermite<Poly>::Hermite() : data({}), m_poly(Poly()), divided_diff({}), x_lists({}), y_lists({}) {}

template <class Poly>
Hermite<Poly>::Hermite(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<int> &nData)
{
    x_lists = xData;
    y_lists = yData;
    data.clear();
    int n = xData.size();
    for (int i = 0; i < n; i++)
    {
        data.push_back(std::make_tuple(xData[i], yData[i], nData[i]));
    }
    // 对data排序，按照xData的大小排序，在此基础上，nData的大小排序
    std::sort(data.begin(), data.end(), [](const std::tuple<double, double, int> &a, const std::tuple<double, double, int> &b) -> bool{
        if (std::abs(std::get<0>(a) - std::get<0>(b)) < 1e-16)
        {
            return std::get<2>(a) < std::get<2>(b);
        }
        return std::get<0>(a) < std::get<0>(b);
    });
    m_poly = interpolate(xData, yData);
}

template <class Poly>
Poly Hermite<Poly>::interpolate(const std::vector<double> &xData, const std::vector<double> &yData)
{
    int n = xData.size();
    divided_diff.resize(n);
    for (int i = 0; i < n; i++)
    {
        divided_diff[i].resize(i + 1, nan(""));
        int index = i-std::get<2>(data[i]);
        divided_diff[i][0] = std::get<1>(data[index]); 
        double divisor = 1.0;
        for (int j = 1; j<=std::get<2>(data[i]); j++)
        {
            divided_diff[i][j] = std::get<1>(data[index+j]) / divisor;
            divisor *= (j+1);
        }
    }
    for (int j = 1; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            if (std::isnan(divided_diff[i][j]))
            {
                divided_diff[i][j] = (divided_diff[i][j - 1] - divided_diff[i - 1][j - 1]) / (std::get<0>(data[i]) - std::get<0>(data[i - j]));
            }
            else
            {
                continue;
            }
        }
    }
    Poly result;
    for (int i = 0; i < n; i++)
    {
        Poly tmp({divided_diff[i][i]});
        for (int j = 0; j < i; j++)
        {
            tmp = tmp * Poly({-std::get<0>(data[j]), 1});
        }
        result = result + tmp;
    }
    return result;
}

template <class Poly>
Poly Hermite<Poly>::add_point(double x, double y, int n)
{
    x_lists.push_back(x);
    y_lists.push_back(y);
    data.push_back(std::make_tuple(x, y, n));
    return interpolate(x_lists, y_lists);
}

template <class Poly>
double Hermite<Poly>::operator()(double x) const
{
    return m_poly(x);
}

template <class Poly>
int Hermite<Poly>::degree() const
{
    return m_poly.get_degree();
}

// template <class P>
// std::ostream &operator<<(std::ostream &os, const Hermite<P> &hermite)
// {
//     os << hermite.m_poly;
//     return os;
// }

template <class Poly>
double Hermite<Poly>::derivative(double x, int n) const
{
    Poly poly = m_poly;
    for (int i = 0; i < n; i++)
    {
        poly = poly.derivative();
    }
    return poly(x);
}

template <class Poly>
double Hermite<Poly>::integral(double a, double b) const
{
    Poly poly = m_poly;
    poly = poly.integral();
    return poly(b) - poly(a);
}

template <class Poly>
Poly Hermite<Poly>::get_polynomial(int n) const
{
    Poly poly = m_poly;
    for (int i = 0; i < n; i++)
    {
        poly = poly.derivative();
    }
    return poly;
}

template class Hermite<Polynomial>;