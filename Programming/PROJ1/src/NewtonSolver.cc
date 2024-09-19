/**
 * @file NewtonSolver.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "NewtonSolver.h"
#include <cmath>

/**
 * @brief Construct a new Newton Solver:: Newton Solver object
 * 
 * @param f The function to be solved
 * @param df The derivative of the function
 * @param x0 The initial value
 */
NewtonSolver::NewtonSolver(double (*f)(double), double (*df)(double), double x0) :
    f(f), df(df), x0(x0), eps(1e-12), MaxIter(50), iter(0), x(0)
{}

/**
 * @brief Construct a new Newton Solver:: Newton Solver object
 * 
 * @param f The function to be solved
 * @param df The derivative of the function
 * @param x0 The initial value
 * @param MaxIter The maximum number of iterations
 */
NewtonSolver::NewtonSolver(double (*f)(double), double (*df)(double), double x0, int MaxIter) :
    f(f), df(df), x0(x0), eps(1e-12), MaxIter(MaxIter), iter(0), x(0)
{}

/**
 * @brief The solve function of the Newton method
 * 
 * The Newton method is a root-finding algorithm that uses the first few terms of the Taylor series of a function \f$f(x)\f$ in the vicinity of a point \f$x_0\f$ to find the root of the function.
 * 
 * The iteration formula is:
 * \f[
 * x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
 * \f]
 * 
 * @return double 
 */
double NewtonSolver::solve()
{
    double fx = f(x0);
    double dfx = df(x0);
    for (iter=0; iter<MaxIter; iter++)
    {
        if (std::abs(fx) < eps)
            break;
        x0 -= fx / dfx; // x0 = x0 - fx / dfx;
        fx = f(x0);
        dfx = df(x0);
    }
    x = x0;
    return x;
}