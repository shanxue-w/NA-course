/**
 * @file Newton.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The implementation of Newton Interpolation.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "Newton.hpp"
#include "Polynomial.hpp"
#include "NewtonPoly.hpp"

/**
 * @brief Construct a new Newton:: Newton object
 */
Newton::Newton()
{
    x_lists = {};
    y_lists = {};
    m_poly = NewtonPoly();
    divided_diff = {};
}

/**
 * @brief Construct a new Newton:: Newton object
 * 
 * @param xData the x values of the data points, distinct, order not required
 * @param yData the y values of the data points
 * 
 * @code {.cc}
 * std::vector<double> xData = {0, 1, 2, 3};
 * std::vector<double> yData = {0, 1, 4, 9};
 * Newton newton(xData, yData); 
 * @endcode
 */
Newton::Newton(const std::vector<double> &xData, const std::vector<double> &yData)
{
    x_lists = xData;
    y_lists = yData;
    m_poly = interpolate(xData, yData);
}

/**
 * @brief Construct a new Newton:: Newton object
 * 
 * @param xData the x values of the data points
 * @param f the function pointer, double
 * 
 * @code {.cc}
 * Newton newton(xData, f);
 * @endcode
 */
Newton::Newton(const std::vector<double> &xData, FuncPtr f)
{
    x_lists = xData;
    y_lists = {};
    for (auto x : xData)
    {
        y_lists.push_back(f(x));
    }
    m_poly = interpolate(xData, y_lists);
}

/**
 * @brief Implement the Newton Interpolation. Based on the divided difference table.
 * 
 * But rarely used outside the class.
 * 
 * @param xData the x values of the data points
 * @param yData the y values of the data points
 * @return NewtonPoly 
 * 
 * @code {.cc}
 * std::vector<double> xData = {0, 1, 2, 3};
 * std::vector<double> yData = {0, 1, 4, 9};
 * NewtonPoly newtonPoly = interpolate(xData, yData);
 * @endcode
 */
NewtonPoly Newton::interpolate(const std::vector<double> &xData, const std::vector<double> &yData)
{
    int n = xData.size();
    divided_diff.resize(n);
    for (int i = 0; i < n; i++)
    {
        divided_diff[i].resize(i + 1);
        divided_diff[i][0] = yData[i];
    }
    for (int j = 1; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            divided_diff[i][j] = (divided_diff[i][j - 1] - divided_diff[i - 1][j - 1]) / (xData[i] - xData[i - j]);
        }
    }
    std::vector<double> w_lists;
    for (int i = 0; i < n; i++)
    {
        w_lists.push_back(divided_diff[i][i]);
    }
    NewtonPoly newtonPoly(xData, yData, w_lists);
    return newtonPoly;
}

/**
 * @brief Add a new point to the Newton Interpolation.
 * 
 * @param x the x value of the new point
 * @param y the y value of the new point
 * @return NewtonPoly 
 * 
 * @code {.cc}
 * NewtonPoly newtonPoly = newton.add_point(1, 2);
 * @endcode
 */
NewtonPoly Newton::add_point(double x, double y)
{
    x_lists.push_back(x);
    y_lists.push_back(y);
    int n = x_lists.size();
    divided_diff.push_back({});
    divided_diff[n - 1].resize(n);
    for (int j = 1; j < n; j++)
    {
        divided_diff[n - 1][j] = (divided_diff[n - 1][j - 1] - divided_diff[n - 2][j - 1]) / (x - x_lists[n - 1 - j]);
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
 * @brief Overload the operator() to get the value of the Newton Interpolation.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * double y = newton(1.0);
 * @endcode
 */
double Newton::operator()(double x) const
{
    return m_poly(x);
}

/**
 * @brief get the degree of the Newton Interpolation.
 * 
 * @return int 
 * 
 * @code {.cc}
 * int n = newton.degree();
 * @endcode
 */
int Newton::degree() const
{
    return m_poly.get_degree();
}

/**
 * @brief derivative of the Newton Interpolation at x.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * double y = newton.derivative(1.0); // p'(1.0)
 * @endcode
 */
double Newton::derivative(double x) const
{
    return m_poly.derivative(x);
}

/**
 * @brief integral of the Newton Interpolation from a to b.
 * 
 * @param a the lower bound of the integral
 * @param b the upper bound of the integral
 * @return double 
 * 
 * @code {.cc}
 * double y = newton.integral(0.0, 1.0); // \f$\int_0^1 p(x)dx\f$
 * @endcode
 */
double Newton::integral(double a, double b) const
{
    return m_poly.integral(a, b);
}

/**
 * @brief get the Newton Polynomial of the Newton Interpolation.
 * 
 * @return NewtonPoly 
 * 
 * @code {.cc}
 * NewtonPoly newtonPoly = newton.getNewtonPoly();
 * @endcode
 */
NewtonPoly Newton::getNewtonPoly() const
{
    return m_poly;
}

/**
 * @brief Convert the Newton Interpolation to Polynomial.
 * 
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly = newton.Convert_to_Polynomial();
 * @endcode
 */
Polynomial Newton::Convert_to_Polynomial() const
{
    return m_poly.Convert_to_Polynomial();
}