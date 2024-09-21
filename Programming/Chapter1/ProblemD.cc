/**
 * @file D.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Qestion D. Testing the SecantSolver.
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "SecantSolver.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

/**
 * @brief for function \f$f(x) = \sin (\frac{x}{2}) - 1\f$
 * 
 * @param x 
 * @return double 
 */
double f1(double x)
{
    return std::sin(x/2.0) - 1.0;
}

/**
 * @brief for function \f$f(x) = e^x - \tan x\f$
 * 
 * @param x 
 * @return double 
 */
double f2(double x)
{
    return std::exp(x) - std::tan(x);
}

/**
 * @brief for function \f$f(x) = x^3-12x^2+3x+1\f$
 * 
 * @param x 
 * @return double 
 */
double f3(double x)
{
    return x*x*x - 12.0*x*x + 3.0*x + 1.0;
}

int main(void)
{
    std::cout << "=========================================\n"
              << "Quesiton D. Testing the SecantSolver.\n\n";
    
    /*<< ================= Part 1 ================= */
    {
        Function F1(f1);
        /**<< for x_0=0.0, x_1=\pi/2 */
        SecantSolver solver1(F1, 0.0, M_PI_2);
        double x1 = solver1.solve();
        std::cout << "The root of f1(x) = sin(x/2) - 1 with x0=0, x1=pi/2 is: " 
                << std::fixed << std::setprecision(9) << x1 
                << "\nf1(x*) = "
                << std::fixed << std::setprecision(9) << f1(x1) << "\n\n";
        
        /**<< for x_0=3\pi, x_1=2\pi */
        SecantSolver solver1_1(F1, 3*M_PI, 4*M_PI);
        double x1_1 = solver1_1.solve();
        std::cout << "The root of f1(x) = sin(x/2) - 1 with x0=3pi, x1=4pi is: " 
                << std::fixed << std::setprecision(9) << x1_1 
                << "\nf1(x*) = "
                << std::fixed << std::setprecision(9) << f1(x1_1) << "\n\n";
    }
    /*<< ================= Part 1 ================= */

    /*<< ================= Part 2 ================= */
    {
        Function F2(f2);
        /**<< for x_0=1.0, x_1=1.4 */
        SecantSolver solver2(F2, 1.0, 1.4);
        double x2 = solver2.solve();
        std::cout << "The root of f2(x) = exp(x) - tan(x) with x0=1.0, x1=1.4 is: " 
                << std::fixed << std::setprecision(9) << x2 
                << "\nf2(x*) = "
                << std::fixed << std::setprecision(9) << f2(x2) << "\n\n";
        
        /**<< for x_0=-4, x_1=-3 */
        SecantSolver solver2_2(F2, -4.0, -3.0);
        double x2_2 = solver2_2.solve();
        std::cout << "The root of f2(x) = exp(x) - tan(x) with x0=-4.0, x1=-3.0 is: " 
                << std::fixed << std::setprecision(9) << x2_2 
                << "\nf2(x*) = "
                << std::fixed << std::setprecision(9) << f2(x2_2) << "\n\n";
    }
    /*<< ================= Part 2 ================= */

    /*<< ================= Part 3 ================= */
    {
        Function F3(f3);
        /**<< for x_0=0.0, x_1=-0.5 */
        SecantSolver solver3(F3, 0.0, -0.5);
        double x3 = solver3.solve();
        std::cout << "The root of f3(x) = x^3 - 12x^2 + 3x + 1 with x0=0.0, x1=-0.5 is: " 
                << std::fixed << std::setprecision(9) << x3 
                << "\nf3(x*) = "
                << std::fixed << std::setprecision(9) << f3(x3) << "\n\n";

        /**<< for x_0=10, x_1=12 */
        SecantSolver solver3_3(F3, 10.0, 12.0);
        double x3_3 = solver3_3.solve();
        std::cout << "The root of f3(x) = x^3 - 12x^2 + 3x + 1 with x0=10.0, x1=12.0 is: " 
                << std::fixed << std::setprecision(9) << x3_3 
                << "\nf3(x*) = "
                << std::fixed << std::setprecision(9) << f3(x3_3) << "\n\n";
    }
    /*<< ================= Part 3 ================= */

    std::cout << "\nQuesiton D. Testing the SecantSolver. Done.\n"
              << "=========================================\n";
    return 0;
}