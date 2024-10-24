/**
 * @file LagPoly.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implementation of Lagrange Polynomial.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "LagPoly.hpp"
#include "NewtonSolver.hpp"
#include <numeric>

/**
 * @brief Construct a new Lag Poly:: Lag Poly object
 */
LagPoly::LagPoly() : x_lists({0}), y_lists({0}), w_lists({0}) {}

/**
 * @brief Construct a new Lag Poly:: Lag Poly object
 * 
 * @param xData 
 * @param yData 
 * 
 * @return None
 * 
 * @code {.cpp}
 * std::vector<double> xData = {1, 2, 3, 4, 5};
 * std::vector<double> yData = {1, 4, 9, 16, 25};
 * LagPoly lagpoly(xData, yData); 
 * @endcode
 * 
 */
LagPoly::LagPoly(const std::vector<double> xData, const std::vector<double> yData) : x_lists(xData), y_lists(yData)
{
    w_lists.resize(xData.size(), 1);
    int n = xData.size();
    for (auto i = 0; i < n; i++)
    {
        for (auto j = 0; j < n; j++)
        {
            if (j != i)
            {
                w_lists[i] /= (xData[i] - xData[j]);
            }
        }
    }
}

/**
 * @brief Construct a new Lag Poly:: Lag Poly object. When we know the weight of each term.
 * 
 * @param xData 
 * @param yData 
 * @param wData 
 * 
 * @code {.cc}
 * std::vector<double> xData = {1, 2, 3, 4, 5};
 * std::vector<double> yData = {1, 4, 9, 16, 25};
 * std::vector<double> wData = {1, 1, 1, 1, 1};
 * LagPoly lagpoly(xData, yData, wData);
 * @endcode
 * 
 */
LagPoly::LagPoly(const std::vector<double> xData, const std::vector<double> yData, const std::vector<double> wData) : 
x_lists(xData), y_lists(yData), w_lists(wData) 
{}

/**
 * @brief Overload the assignment operator.
 * 
 * @param rhs 
 * @return LagPoly& 
 * 
 * @code {.cc}
 * LagPoly lagpoly1;
 * LagPoly lagpoly2;
 * lagpoly1 = lagpoly2;
 * @endcode
 * 
 */
LagPoly &LagPoly::operator=(const LagPoly &rhs)
{
    if (this != &rhs)
    {
        x_lists = rhs.x_lists;
        y_lists = rhs.y_lists;
        w_lists = rhs.w_lists;
    }
    return *this;
}

/**
 * @brief Overload the operator() to get the value of the Lagrange Polynomial.
 * 
 * @param x 
 * @return double 
 * 
 * @code {.cc}
 * LagPoly lagpoly;
 * double y = lagpoly(1.0);
 * @endcode 
 */
double LagPoly::operator()(double x) const
{
    // 分子 分母的英文 为 numerator 和 denominator
    std::vector<double> numerator(x_lists.size(), 0);
    std::vector<double> denominator(x_lists.size(), 0);
    int n = x_lists.size();
    for (auto i=0; i<n; i++)
    {
        double tmp = w_lists[i] / (x - x_lists[i]);
        numerator[i] = tmp * y_lists[i];
        denominator[i] = tmp;
    }
    return std::accumulate(numerator.begin(), numerator.end(), 0.0) / std::accumulate(denominator.begin(), denominator.end(), 0.0);
}

/**
 * @brief The derivative of the Lagrange Polynomial.
 * 
 * @param x 
 * @return double 
 * 
 * @code {.cc}
 * LagPoly lagpoly;
 * double y = lagpoly.derivative(1.0);
 * @endcode
 * 
 */
double LagPoly::derivative(double x) const
{
    double epison = 1e-4;
    return (this->operator()(x+epison) - this->operator()(x-epison)) / (2*epison);
}

/**
 * @brief The integral of the Lagrange Polynomial.
 * 
 * @param a the lower bound
 * @param b the upper bound
 * @return double 
 * 
 * @code {.cc}
 * LagPoly lagpoly;
 * double y = lagpoly.integral(0.0, 1.0);
 * @endcode
 * 
 */
double LagPoly::integral(double a, double b) const
{
    double delta_x = 1e-3;
    double result = 0.0;
    for (double x = a; x < b; x += delta_x)
    {
        result += this->operator()(x) * delta_x;
    }
    return result;
}

/**
 * @brief The degree of the Lagrange Polynomial.
 * 
 * @return int 
 * 
 * @code {.cc}
 * LagPoly lagpoly;
 * int n = lagpoly.get_degree();
 * @endcode
 */
int LagPoly::get_degree() const
{
    return x_lists.size() - 1;
}

/**
 * @brief Overload the operator<< to output the Lagrange Polynomial.
 * 
 * @param os 
 * @param lagpoly 
 * @return std::ostream& 
 * 
 * @code {.cc}
 * LagPoly lagpoly;
 * std::cout << lagpoly << std::endl;
 * @endcode
 */
std::ostream &operator<<(std::ostream &os, const LagPoly &lagpoly)
{
    os << "Lagrange Polynomial: ";
    int n = lagpoly.x_lists.size();
    os << "[";
    for (auto i = 0; i<n; i++)
    {
        if (i==0)
        {
            os << lagpoly.y_lists[i];
        }
        else
        {
            if (lagpoly.y_lists[i] >= 0)
            {
                os << " + " << lagpoly.y_lists[i];
            }
            else
            {
                os << " - " << -lagpoly.y_lists[i];
            }
        }
        os << " * (" << lagpoly.w_lists[i] << ") / (x";
        if (lagpoly.x_lists[i] >= 0)
        {
            os << " - " << lagpoly.x_lists[i] << ")";
        }
        else
        {
            os << " + " << -lagpoly.x_lists[i] << ")";
        }
    }
    os << "] / [";
    for (auto i = 0; i<n; i++)
    {
        if (i == 0)
        {
            os << "( " << lagpoly.w_lists[i] << ") / (x";
            if (lagpoly.x_lists[i] >= 0)
            {
                os << " - " << lagpoly.x_lists[i] << ")";
            }
            else
            {
                os << " + " << -lagpoly.x_lists[i] << ")";
            }
        }
        else
        {
            os << "+( " << lagpoly.w_lists[i] << ") / (x";
            if (lagpoly.x_lists[i] >= 0)
            {
                os << " - " << lagpoly.x_lists[i] << ")";
            }
            else
            {
                os << " + " << -lagpoly.x_lists[i] << ")";
            }
        }
    }
    os << "]";
    return os;
}