#include "Hermite.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include "Polynomial.hpp"
#include "NewtonPoly.hpp"

/**
 * @brief Construct a new Hermite:: Hermite object
 */
Hermite::Hermite() : data({}), divided_diff({}), x_lists({}), y_lists({}) {}

/**
 * @brief Construct a new Hermite:: Hermite object
 * 
 * @param xData the x values of the data points
 * @param yData the y values of the data points, associated with nData, means the n-th derivative of f(x) at xData[i]
 * @param nData represent the n-th derivative of f(x) at xData[i]
 * 
 * @code {.cc}
 * std::vector<double> xData = {0, 1, 0, 1};
 * std::vector<double> yData = {0, 1, 3, 4};
 * std::vector<int> nData = {0, 0, 1, 1};
 * // means f(0) = 0, f(1) = 1, f'(0) = 3, f'(1) = 4
 * Hermite hermite(xData, yData, nData);
 * @endcode
 */
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

/**
 * @brief Implement the Hermite Interpolation. Based on the divided difference table.
 * 
 * But rarely used outside the class.
 * 
 * @param data_lists the data lists, each element is a tuple (x, y, n), n means the n-th derivative of f(x) at x
 * @return NewtonPoly 
 * 
 * @code {.cc}
 * std::vector<std::tuple<double, double, int>> data_lists = {std::make_tuple(0, 0, 0), std::make_tuple(1, 1, 0), std::make_tuple(0, 3, 1), std::make_tuple(1, 4, 1)};
 * Hermite hermite;
 * NewtonPoly newtonPoly = hermite.interpolate(data_lists);
 * @endcode
 */
NewtonPoly Hermite::interpolate(std::vector<std::tuple<double, double, int>> &data_lists)
{
    int n = x_lists.size();
    divided_diff.resize(n);
    for (int i = 0; i < n; i++)
    {
        divided_diff[i].resize(i + 1, nan("")); /**set the unknown to nan */
        int index = i-std::get<2>(data_lists[i]);
        divided_diff[i][0] = std::get<1>(data_lists[index]); 
        double divisor = 1.0;
        for (int j = 1; j<=std::get<2>(data_lists[i]); j++)
        {
            divided_diff[i][j] = std::get<1>(data_lists[index+j]) / divisor;
            divisor *= (j+1);
        }
    }
    /** Build the difference table */
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

/**
 * @brief Add a new point in the Hermite Interpolation.
 * 
 * @param x the x value of the new point
 * @param y the y value of the new point, associated with n, means the n-th derivative of f(x) at x
 * @param n represent the n-th derivative of f(x) at x
 * @return NewtonPoly 
 * 
 * @code {.cc}
 * NewtonPoly newtonPoly = hermite.add_point(1, 2, 0);
 * @endcode
 */
NewtonPoly Hermite::add_point(double x, double y, int n)
{
    x_lists.push_back(x);
    y_lists.push_back(y);
    data.push_back(std::make_tuple(x, y, n));
    NewtonPoly newtonPoly;
    newtonPoly = interpolate(data);
    return newtonPoly;
}

/**
 * @brief Overload the operator() to get the value of the Hermite Interpolation.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * double y = hermite(1.0);
 * @endcode
 */
double Hermite::operator()(double x) const
{
    return m_poly(x);
}

/**
 * @brief get the degree of the Hermite Interpolation.
 * 
 * @return int 
 * @code {.cc}
 * int n = hermite.degree();
 * @endcode
 */
int Hermite::degree() const
{
    return m_poly.get_degree();
}


/**
 * @brief derivative of the Hermite Interpolation at x.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * double y = hermite.derivative(1.0); // p'(1.0)
 * @endcode
 */
double Hermite::derivative(double x) const
{
    return m_poly.derivative(x);
}

/**
 * @brief integral of the Hermite Interpolation from a to b.
 * 
 * @param a the lower bound of the integral
 * @param b the upper bound of the integral
 * @return double 
 * 
 * @code {.cc}
 * double y = hermite.integral(0.0, 1.0); // \f$\int_0^1 p(x)dx\f$
 * @endcode
 */
double Hermite::integral(double a, double b) const
{
    return m_poly.integral(a, b);
}

/**
 * @brief get the Newton Polynomial of the Hermite Interpolation.
 * 
 * @return NewtonPoly 
 * 
 * @code {.cc}
 * NewtonPoly newtonPoly = hermite.get_polynomial();
 * @endcode
 */
NewtonPoly Hermite::get_polynomial()
{
    return m_poly;
}

/**
 * @brief convert the Hermite Interpolation to Polynomial.
 * 
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly = hermite.Convert_to_Polynomial();
 * @endcode
 */
Polynomial Hermite::Convert_to_Polynomial() const
{
    return m_poly.Convert_to_Polynomial();
}