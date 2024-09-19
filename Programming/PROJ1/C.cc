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

#include "NewtonSolver.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

double f(double x)
{
    return x - std::tan(x);
}

double df(double x)
{
    double cosx = std::cos(x);
    return 1.0 - 1.0/(cosx * cosx);
}

int main(void)
{
    NewtonSolver solver1(f, df, 4.5);
    double x1 = solver1.solve();
    std::cout << "The root of f(x) = x - tan(x) near 4.5 is: " 
            << std::fixed << std::setprecision(9) << x1 
            << "\nf(x1) = "
            << std::fixed << std::setprecision(9) << f(x1) << "\n\n";
    
    NewtonSolver solver2(f, df, 7.7);
    double x2 = solver2.solve();
    std::cout << "The root of f(x) = x - tan(x) near 7.7 is: " 
            << std::fixed << std::setprecision(9) << x2 
            << "\nf(x2) = "
            << std::fixed << std::setprecision(9) << f(x2) << "\n\n";

    return 0;
}