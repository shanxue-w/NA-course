/**
 * @file F.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Quesiton F
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

double l;
double beta_1;
double h;
double D;

double A, B, C, E;
const double M_OVER_180 = M_PI/180.0; // 180/M_PI

double f(double x)
{
    x = x * M_OVER_180;
    double sinx = std::sin(x);
    double cosx = std::cos(x);
    return A*sinx*cosx + B*sinx*sinx - C*cosx - E*sinx; 
}

double df(double x)
{
    x = x * M_OVER_180;
    double sinx = std::sin(x);
    double cosx = std::cos(x);
    return M_OVER_180* (A*cosx*cosx - A*sinx*sinx + 2*B*sinx*cosx + C*sinx - E*cosx);
}

int main(void)
{
    // Part 1
    {
        l = 89;
        beta_1 = 11.5/180.0*M_PI;
        h = 49;
        D = 55;
        A = l*std::sin(beta_1);
        B = l*std::cos(beta_1);
        C = (h+0.5*D)*std::sin(beta_1) - 0.5*D*std::tan(beta_1);
        E = (h+0.5*D)*std::cos(beta_1) - 0.5*D;
        NewtonSolver solver1(f, df, 20);
        double x1 = solver1.solve();
        std::cout << "When l=89, h=49, beta_1=11.5, D=55, the angle is: " 
                << std::fixed << std::setprecision(9) << x1 
                << "\nf(x1) = "
                << std::fixed << std::setprecision(9) << f(x1) << "\n\n";
    }


    // Part 2
    {
        l = 89;
        beta_1 = 11.5/180.0*M_PI;
        h = 49;
        D = 30;
        A = l*std::sin(beta_1);
        B = l*std::cos(beta_1);
        C = (h+0.5*D)*std::sin(beta_1) - 0.5*D*std::tan(beta_1);
        E = (h+0.5*D)*std::cos(beta_1) - 0.5*D;

        NewtonSolver solver2(f, df, 33);
        double x2 = solver2.solve();
        std::cout << "When l=89, h=49, beta_1=11.5, D=30, the angle is: " 
                << std::fixed << std::setprecision(9) << x2 
                << "\nf(x2) = "
                << std::fixed << std::setprecision(9) << f(x2) << "\n\n";
    }

    // Part 3
    {
        l = 89;
        beta_1 = 11.5/180.0*M_PI;
        h = 49;
        D = 30;
        A = l*std::sin(beta_1);
        B = l*std::cos(beta_1);
        C = (h+0.5*D)*std::sin(beta_1) - 0.5*D*std::tan(beta_1);
        E = (h+0.5*D)*std::cos(beta_1) - 0.5*D;

        SecantSolver solver3(f, 0, 33);
        double x3 = solver3.solve();
        std::cout << "When l=89, h=49, beta_1=11.5, D=30, the angle is: " 
                << std::fixed << std::setprecision(9) << x3 
                << "\nf(x3) = "
                << std::fixed << std::setprecision(9) << f(x3) << "\n\n";
    }
    
    return 0;
}