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

#include "SecantSolver.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

double f1(double x)
{
    return std::sin(x/2.0) - 1.0;
}

double f2(double x)
{
    return std::exp(x) - std::tan(x);
}

double f3(double x)
{
    return x*x*x - 12.0*x*x + 3.0*x + 1.0;
}

int main(void)
{
    SecantSolver solver1(f1, 0.0, M_PI_2);
    double x1 = solver1.solve();
    std::cout << "The root of f1(x) = sin(x/2) - 1 with x0=0, x1=pi/2 is: " 
            << std::fixed << std::setprecision(9) << x1 
            << "\nf1(x1) = "
            << std::fixed << std::setprecision(9) << f1(x1) << "\n\n";
    
    SecantSolver solver2(f2, 1.0, 1.4);
    double x2 = solver2.solve();
    std::cout << "The root of f2(x) = exp(x) - tan(x) with x0=1.0, x1=1.4 is: " 
            << std::fixed << std::setprecision(9) << x2 
            << "\nf2(x2) = "
            << std::fixed << std::setprecision(9) << f2(x2) << "\n\n";
    
    SecantSolver solver3(f3, 0.0, -0.5);
    double x3 = solver3.solve();
    std::cout << "The root of f3(x) = x^3 - 12x^2 + 3x + 1 with x0=0.0, x1=-0.5 is: " 
            << std::fixed << std::setprecision(9) << x3 
            << "\nf3(x3) = "
            << std::fixed << std::setprecision(9) << f3(x3) << "\n\n";

    return 0;
}