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

#include "BisectionSolver.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>

/**
 * @brief Construct a new Bisection Solver:: Bisection Solver object
 * 
 * @param f The function to be solved
 * @param a The left end of the interval
 * @param b The right end of the interval
 */
BisectionSolver::BisectionSolver(double (*f)(double), double a, double b) :
    f(f), a(a), b(b), eps(1e-12), MaxIter(100), iter(0), x(0)
{}

/**
 * @brief Construct a new Bisection Solver:: Bisection Solver object
 * 
 * @param f The function to be solved
 * @param a The left end of the interval
 * @param b The right end of the interval
 * @param MaxIter The maximum number of iterations
 */
BisectionSolver::BisectionSolver(double (*f)(double), double a, double b, int MaxIter) :
    f(f), a(a), b(b), eps(1e-12), MaxIter(MaxIter), iter(0), x(0)
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
        if (std::abs(fc) < eps)
            break;
        else if (fa * fc > 0)
            a = c;
    }
    x = c;
    return x;
}