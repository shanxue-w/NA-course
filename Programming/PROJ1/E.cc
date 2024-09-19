/**
 * @file E.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Question E
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "BisectionSolver.h"
#include "NewtonSolver.h"
#include "SecantSolver.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

double f(double h)
{
    return 10*(M_PI_2 - std::asin(h) - h*std::sqrt(1.0 - h*h)) - 12.4;
}

double df(double h)
{
    return 10*(-1.0/std::sqrt(1.0 - h*h) - std::sqrt(1.0 - h*h) + h*h/std::sqrt(1.0 - h*h));
}

int main(void)
{
    BisectionSolver solver1(f, 0.0, 1.0);
    double h1 = solver1.solve();
    if (solver1.getIter() >= 0)
    {
        std::cout << "The root of f(h) = 10*(pi/2 - asin(h) - h*sqrt(1-h^2)) - 12.4 is: " 
                  << std::fixed << std::setprecision(9) << h1 
                  << "\nf(h1) = "
                  << std::fixed << std::setprecision(9) << f(h1) << "\n\n";
    }

    NewtonSolver solver2(f, df, 0.5);
    double h2 = solver2.solve();
    std::cout << "The root of f(h) = 10*(pi/2 - asin(h) - h*sqrt(1-h^2)) - 12.4 is: " 
              << std::fixed << std::setprecision(9) << h2 
              << "\nf(h2) = "
              << std::fixed << std::setprecision(9) << f(h2) << "\n\n";

    SecantSolver solver3(f, 0.0, 1.0);
    double h3 = solver3.solve();
    std::cout << "The root of f(h) = 10*(pi/2 - asin(h) - h*sqrt(1-h^2)) - 12.4 is: " 
              << std::fixed << std::setprecision(9) << h3 
              << "\nf(h3) = "
              << std::fixed << std::setprecision(9) << f(h3) << "\n\n";
    return 0;
}