/**
 * @file Lagrange.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The implementation of Lagrange Interpolation.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "LagPoly.hpp"
#include "Lagrange.hpp"

/**
 * @brief Construct a new Lagrange:: Lagrange object
 * 
 */
Lagrange::Lagrange() : poly() {}

/**
 * @brief Construct a new Lagrange:: Lagrange object
 * 
 * @param xData the x values of the data points
 * @param yData the y values of the data points
 * 
 * @code {.cc}
 * std::vector<double> xData = {0, 1, 2, 3};
 * std::vector<double> yData = {0, 1, 4, 9};
 * Lagrange lagrange(xData, yData);
 * @endcode
 */
Lagrange::Lagrange(const std::vector<double> &xData, const std::vector<double> &yData) : poly(xData, yData) {}

/**
 * @brief Construct a new Lagrange:: Lagrange object
 * 
 * Used in copy, rarely used outside the class.
 * 
 * @param xData the x values of the data points
 * @param yData the y values of the data points
 * @param wData the weight values of the data points
 * 
 * @code {.cc}
 * std::vector<double> xData = {0, 1, 2, 3};
 * std::vector<double> yData = {0, 1, 4, 9};
 * std::vector<double> wData = {1, 1, 1, 1};
 * Lagrange lagrange(xData, yData, wData);
 * @endcode
 */
Lagrange::Lagrange(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &wData) : 
poly(xData, yData, wData) 
{}

/**
 * @brief Overload the operator() to get the value of the Lagrange Interpolation.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * double y = lagrange(1.0);
 * @endcode
 */
double Lagrange::operator()(double x) const
{
    return poly(x);
}

/**
 * @brief get the degree of the Lagrange Interpolation.
 * 
 * @return int 
 * @code {.cc}
 * int n = lagrange.degree();
 * @endcode
 */
int Lagrange::degree() const
{
    return poly.get_degree();
}

/**
 * @brief get the polynomial of the Lagrange Interpolation.
 * 
 * @return LagPoly 
 * @code {.cc}
 * LagPoly poly = lagrange.get_polynomial();
 * @endcode
 */
LagPoly Lagrange::get_polynomial() const
{
    return poly;
}

/**
 * @brief Overload the operator<< to output the Lagrange Interpolation.
 * 
 * @param os the output stream
 * @param lagrange the Lagrange object
 * @return std::ostream& 
 * 
 * @code {.cc}
 * std::cout << lagrange << std::endl;
 * @endcode
 */
std::ostream &operator<<(std::ostream &os, const Lagrange &lagrange)
{
    os << lagrange.poly;
    return os;
}
