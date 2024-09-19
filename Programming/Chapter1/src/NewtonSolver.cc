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

#include "NewtonSolver.hpp"


/**
 * @brief Construct a new Newton Solver:: Newton Solver object
 * 
 * @param f The function to be solved and its derivative
 * @param x0 The initial value
 * @param eps The tolerance of f(x)
 * @param delta The tolerance of the interval
 * @param MaxIter The maximum number of iterations
 */
NewtonSolver::NewtonSolver(Function f, double x0, double eps, double delta, int MaxIter)
    : f(f), x0(x0), eps(eps), delta(delta), MaxIter(MaxIter), iter(0), x(0.0)
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
    double dfx = f.derivative(x0);
    double delta_x = 0;
    for (iter=0; iter<MaxIter; iter++)
    {
        if (std::abs(fx) < eps)
            break;
        delta_x = fx / dfx;
        if (std::abs(delta_x) < delta)
            break;
        x0 -= delta_x; /**<< x0 = x0 - f(x)/df(x) */
        fx = f(x0);
        dfx = f.derivative(x0);
    }
    x = x0;
    return x;
}