/**
 * @file BisectionSolver.c
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "BisectionSolver.hpp"

/**
 * @brief Construct a new Bisection Solver:: Bisection Solver object
 * 
 * @param f The function to be solved
 * @param a The left end of the interval
 * @param b The right end of the interval
 * @param eps The tolerance of f(x)
 * @param delta The tolerance of the interval
 * @param MaxIter The maximum number of iterations
 */
BisectionSolver::BisectionSolver(Function f, double a, double b, double eps, double delta, int MaxIter) :
    f(f), a(a), b(b), eps(eps), delta(delta), MaxIter(MaxIter), iter(0), x(0)
{}


/**
 * @brief The solve function of the Bisection method
 * 
 * @return double 
 */
double BisectionSolver::solve()
{
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0)
    {
        std::cerr << "The interval does not satisfy the intermediate value theorem! Please check the interval." << std::endl;
        iter = -1;
        x = 0;
        return x;
    }
    if (a > b)
        std::swap(a, b);
    double h = b - a;
    double c=0, fc=0;
    for (iter=0; iter<MaxIter; iter++)
    {
        h /= 2;
        c = a + h;
        fc = f(c);
        if (std::abs(fc) < eps || h < delta)
            break;
        else if (fa * fc > 0)
            a = c;
    }
    if (std::abs(fc) >= 1e2)
    {
        std::cerr << "The iteration does not converge! Please check the interval. Or may be the function is not continuous." << std::endl;
        iter = -1;
        x = 0;
        return x;
    }
    x = c;
    return x;
}