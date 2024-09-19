/**
 * @file B.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Question B. Testing the BisectionSolver.
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "BisectionSolver.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

/**
 * @brief for function \f$f(x) = \frac{1}{x} - \tan(x)\f$
 * 
 * @param x 
 * @return double 
 */
double f1(double x)
{
    return 1.0/x - std::tan(x);
}

/**
 * @brief for function \f$f(x) = \frac{1}{x} - 2^x\f$
 * 
 * @param x 
 * @return double 
 */
double f2(double x)
{
    return 1.0/x - std::pow(2.0, x);
}

/**
 * @brief for function \f$f(x) = 0.5^x + e^x + 2\cos(x) - 6\f$
 * 
 * @param x 
 * @return double 
 */
double f3(double x)
{
    return std::pow(0.5, x) + std::exp(x) + 2.0*std::cos(x) - 6.0;
}


/**
 * @brief The floating-point error of the Horner algorithm is too large, so it is not used here.
 * for function \f$ f(x) = \frac{x^3 + 4x^2 + 3x + 5}{2x^3 - 9x^2 + 18x - 2} \f$
 * 
 * But to find the root, we don't need to calculate the denominator, so it's deduced to:
 * 
 * \f[
 * f(x) = x^3 + 4x^2 + 3x + 5 = 0
 * \f]
 * 
 * @param x 
 * @return double 
 */
double f4(double x)
{
    return x*x*x + 4.0*x*x + 3*x + 5;
}

int main(void)
{
    BisectionSolver solver1(f1, 0.0, M_PI_2);
    double x1 = solver1.solve();
    if (solver1.getIter() >= 0)
    {
        std::cout << "The root of f1(x) = 1/x - tan(x) is: " 
                << std::fixed << std::setprecision(9) << x1 
                << "\nf1(x1) = "
                << std::fixed << std::setprecision(9) << f1(x1) << "\n\n";
    }

    BisectionSolver solver2(f2, 0.0, 1.0);
    double x2 = solver2.solve();
    if (solver2.getIter() >= 0)
    {
        std::cout << "The root of f2(x) = 1/x - 2^x is: " 
                << std::fixed << std::setprecision(9) << x2
                << "\nf2(x2) = "
                << std::fixed << std::setprecision(9) << f2(x2) << "\n\n";
    }

    BisectionSolver solver3(f3, 1.0, 3.0);
    double x3 = solver3.solve();
    if (solver3.getIter() >= 0)
    {
        std::cout << "The root of f3(x) = 0.5^x + e^x + 2cos(x) - 6 is: "
                << std::fixed << std::setprecision(9) << x3
                << "\nf3(x3) = "
                << std::fixed << std::setprecision(9) << f3(x3) << "\n\n";
    }

    BisectionSolver solver4(f4, 0.0, 4.0);
    double x4 = solver4.solve();
    if (solver4.getIter() >= 0)
    {
        std::cout << "The root of f4(x) = (x^3 + 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2) is: " 
                << std::fixed << std::setprecision(9) << x4
                << "\nf4(x4) = "
                << std::fixed << std::setprecision(9) << f4(x4) << "\n\n";
    }

    return 0;
}