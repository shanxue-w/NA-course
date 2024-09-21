/**
 * @file C.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Question C. Testing the NewtonSolver.
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "NewtonSolver.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

/**
 * @brief Function of this question
 * 
 * \f[
 * f(x) = x - \tan x
 * \f]
 * 
 * @param x 
 * @return double 
 */
double f(double x)
{
    return x - std::tan(x);
}

/**
 * @brief The derivative of the function f.
 * 
 * \f[
 * f^\prime (x) = 1 - \frac{1}{\cos^2 x}
 * \f]
 * 
 * @param x 
 * @return double 
 */
double df(double x)
{
    double cosx = std::cos(x);
    return 1.0 - 1.0/(cosx * cosx);
}

int main(void)
{
    std::cout << "=========================================\n"
              << "Quesiton C. Testing the NewtonSolver.\n\n";
    
    Function F(f, df);

    /*<< ================= Part 1, for root near 4.5 ================= */
    {
        NewtonSolver solver1(F, 4.5);
        double x1 = solver1.solve();
        std::cout << "The root of f(x) = x - tan(x) near 4.5 is: " 
                << std::fixed << std::setprecision(9) << x1 
                << "\nf(x1) = "
                << std::fixed << std::setprecision(9) << f(x1) << "\n\n";
    }
    /*<< ================= Part 1, for root near 4.5 ================= */

    /*<< ================= Part 2, for root near 7.7 ================= */
    {
        NewtonSolver solver2(F, 7.7);
        double x2 = solver2.solve();
        std::cout << "The root of f(x) = x - tan(x) near 7.7 is: " 
                << std::fixed << std::setprecision(9) << x2 
                << "\nf(x2) = "
                << std::fixed << std::setprecision(9) << f(x2) << "\n\n";
    }
    /*<< ================= Part 2, for root near 7.7 ================= */

    std::cout << "\nQuesiton C. Testing the NewtonSolver. Done.\n"
              << "=========================================\n";
    return 0;
}