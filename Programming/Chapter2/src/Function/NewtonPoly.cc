/**
 * @file NewtonPoly.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The class for Newton Polynomial, providing the basic operations of Newton Polynomial.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "NewtonPoly.hpp"
#include "NewtonSolver.hpp"
#include <numeric>

/**
 * @brief Construct a new Newton Poly:: Newton Poly object
 */
NewtonPoly::NewtonPoly() : x_lists({0}), y_lists({0}), w_lists({0}) {}

/**
 * @brief Construct a new Newton Poly:: Newton Poly object
 * 
 * This is mainly used to copy, not for interpolation.
 * 
 * @param xData the x values of the data points
 * @param yData the y values of the data points
 * @param wData the weight of each term
 * 
 * @return None
 * 
 * @code {.cc}
 * std::vector<double> xData = {0, 1, 2};
 * std::vector<double> yData = {0, 1, 4};
 * std::vector<double> wData = {0, 1, 1};
 * NewtonPoly newtonpoly(xData, yData, wData); 
 * @endcode
 * 
 */
NewtonPoly::NewtonPoly(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &wData) : 
x_lists(xData), y_lists(yData), w_lists(wData) 
{}

/**
 * @brief Overload the copy constructor.
 * 
 * @param rhs 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly1;
 * NewtonPoly newtonpoly2(newtonpoly1);
 * @endcode
 * 
 */
NewtonPoly::NewtonPoly(const NewtonPoly &rhs) : x_lists(rhs.x_lists), y_lists(rhs.y_lists), w_lists(rhs.w_lists) {}

/**
 * @brief Overload the assignment operator.
 * 
 * @param rhs the right hand side of the assignment operator
 * @return NewtonPoly& 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly1;
 * NewtonPoly newtonpoly2;
 * newtonpoly1 = newtonpoly2;
 * @endcode
 */
NewtonPoly &NewtonPoly::operator=(const NewtonPoly &rhs)
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
 * @brief Overload the operator() to get the value of the Newton Polynomial.
 * 
 * Use the Horner's method to calculate the value of the Newton Polynomial.
 * \f[
 * f(x) = \Sigma_{i=0}^{n} w_i \Pi_{j=0}^{i-1} (x - x_j) = w_0 + (x-x_0)(w_1 + (x-x_1)(w_2 + \cdots + (x-x_{n-1})w_n) \cdots )
 * \f]
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly;
 * double y = newtonpoly(1.0);
 * @endcode
 */
double NewtonPoly::operator()(double x) const
{
    double result = 0.0;
    int n = x_lists.size();
    for (auto i= n-1; i >= 0; i--)
    {
        result = result * (x - x_lists[i]) + w_lists[i];
    }
    return result;
}

/**
 * @brief Calculate the derivative of the Newton Polynomial.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly;
 * double y = newtonpoly.derivative(1.0);
 * @endcode
 */
double NewtonPoly::derivative(double x) const
{
    double epison = 1e-4;
    return (this->operator()(x+epison) - this->operator()(x-epison)) / (2*epison);
}

/**
 * @brief Calculate the integral of the Newton Polynomial.
 * 
 * @param a the lower bound of the integral
 * @param b the upper bound of the integral
 * @return double 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly;
 * double y = newtonpoly.integral(0.0, 1.0);
 * @endcode
 */
double NewtonPoly::integral(double a, double b) const
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
 * @brief Get the degree of the Newton Polynomial.
 * 
 * @return int 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly;
 * int n = newtonpoly.get_degree();
 * @endcode
 */
int NewtonPoly::get_degree() const
{
    return x_lists.size()-1;
}

/**
 * @brief Overload the operator<< to output the Newton Polynomial.
 * 
 * @param os 
 * @param newtonpoly 
 * @return std::ostream& 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly;
 * std::cout << newtonpoly << std::endl;
 * @endcode
 */
std::ostream &operator<<(std::ostream &os, const NewtonPoly &newtonpoly)
{
    int n = newtonpoly.x_lists.size();
    for (auto i=0; i<n; i++)
    {
        if (i == 0)
        {
            os << newtonpoly.w_lists[i];
        }
        else
        {
            if (newtonpoly.w_lists[i] > 0)
            {
                os << " + " << newtonpoly.w_lists[i];
            }
            else
            {
                os << " - " << -newtonpoly.w_lists[i];
            }
        }
        for (auto j=0; j<i; j++)
        {
            if (newtonpoly.x_lists[j] < 0)
            {
                os << "(x + " << -newtonpoly.x_lists[j] << ")";
            }
            else
            {
                os << "(x - " << newtonpoly.x_lists[j] << ")";
            }
        }
    }
    return os;
}

/**
 * @brief Convert the Newton Polynomial to Polynomial.
 * 
 * @return Polynomial 
 * 
 * @code {.cc}
 * NewtonPoly newtonpoly;
 * Polynomial poly = newtonpoly.Convert_to_Polynomial();
 * @endcode
 */
Polynomial NewtonPoly::Convert_to_Polynomial() const
{
    int n = x_lists.size();
    Polynomial poly({w_lists[0]});
    for (auto i=1; i<n; i++)
    {
        Polynomial temp({w_lists[i]});
        for (auto j=0; j<i; j++)
        {
            Polynomial factor({-x_lists[j], 1});
            temp = temp * factor;
        }
        poly = poly + temp;
    }
    return poly;
}