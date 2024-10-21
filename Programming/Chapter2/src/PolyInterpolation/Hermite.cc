#include "Hermite.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include "Polynomial.hpp"
#include "NewtonPoly.hpp"


Hermite::Hermite() : data({}), divided_diff({}), x_lists({}), y_lists({}) {}


Hermite::Hermite(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<int> &nData)
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
    // std::cout << data << std::endl;
    for (auto i=0; i<n; i++)
    {
        x_lists[i] = std::get<0>(data[i]);
    }
    m_poly = interpolate(data);
}


NewtonPoly Hermite::interpolate(std::vector<std::tuple<double, double, int>> &data_lists)
{
    int n = x_lists.size();
    divided_diff.resize(n);
    for (int i = 0; i < n; i++)
    {
        divided_diff[i].resize(i + 1, nan(""));
        int index = i-std::get<2>(data_lists[i]);
        divided_diff[i][0] = std::get<1>(data_lists[index]); 
        double divisor = 1.0;
        for (int j = 1; j<=std::get<2>(data_lists[i]); j++)
        {
            divided_diff[i][j] = std::get<1>(data_lists[index+j]) / divisor;
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
    std::vector<double> w_lists;
    for (int i = 0; i < n; i++)
    {
        w_lists.push_back(divided_diff[i][i]);
    }
    NewtonPoly newtonPoly(x_lists, y_lists, w_lists);
    return newtonPoly;
}


NewtonPoly Hermite::add_point(double x, double y, int n)
{
    x_lists.push_back(x);
    y_lists.push_back(y);
    data.push_back(std::make_tuple(x, y, n));
    NewtonPoly newtonPoly;
    newtonPoly = interpolate(data);
    return newtonPoly;
}


double Hermite::operator()(double x) const
{
    return m_poly(x);
}


int Hermite::degree() const
{
    return m_poly.get_degree();
}

// template <class P>
// std::ostream &operator<<(std::ostream &os, const Hermite<P> &hermite)
// {
//     os << hermite.m_poly;
//     return os;
// }


double Hermite::derivative(double x) const
{
    return m_poly.derivative(x);
}

double Hermite::integral(double a, double b) const
{
    return m_poly.integral(a, b);
}


NewtonPoly Hermite::get_polynomial()
{
    return m_poly;
}

Polynomial Hermite::Convert_to_Polynomial() const
{
    return m_poly.Convert_to_Polynomial();
}