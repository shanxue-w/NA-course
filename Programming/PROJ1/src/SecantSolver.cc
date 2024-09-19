/**
 * @file SecantSolver.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "SecantSolver.h"
#include <cmath>

/**
 * @brief Construct a new Secant Solver:: Secant Solver object
 * 
 * @param f The function to be solved
 * @param x0 The initial value
 * @param x1 The initial value
 */
SecantSolver::SecantSolver(double (*f)(double), double x0, double x1) :
    f(f), x0(x0), x1(x1), eps(1e-12), MaxIter(50), iter(0), x(0)
{}

/**
 * @brief Construct a new Secant Solver:: Secant Solver object
 * 
 * @param f The function to be solved
 * @param x0 The initial value
 * @param x1 The initial value
 * @param MaxIter The maximum number of iterations
 */
SecantSolver::SecantSolver(double (*f)(double), double x0, double x1, int MaxIter) :
    f(f), x0(x0), x1(x1), eps(1e-12), MaxIter(MaxIter), iter(0), x(0)
{}

/**
 * @brief The solve function of the Secant method
 * 
 * The Secant method is a root-finding algorithm that uses the secant line between two points \f$(x_0, f(x_0))\f$ and \f$(x_1, f(x_1))\f$ to find the root of the function.
 * 
 * The iteration formula is:
 * \f[
 * x_{n+1} = x_n - f(x_n) \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}
 * \f]
 * 
 * @return double 
 */
double SecantSolver::solve()
{
    double fx0 = f(x0);
    double fx1 = f(x1);
    double x2=0;
    for (iter=0; iter<MaxIter; iter++)
    {
        if (std::abs(fx1) < eps)
            break;
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        x0 = x1;
        fx0 = fx1;
        x1 = x2;
        fx1 = f(x1);
    }
    x = x2;
    return x;
}