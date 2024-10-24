/**
 * @file NewtonSolver.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The implementation of Newton Solver. Template class.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "NewtonSolver.hpp"
#include "SolverFunction.hpp"
#include "Polynomial.hpp"

/**
 * @brief Construct a new Newton Solver:: Newton Solver object
 * 
 * Usually use the Polynomial or SolverFunction as the template parameter.
 * 
 * @tparam F the function class
 * @param f the function
 * @param x0 the initial value
 * @param eps the tolerance of the function value
 * @param delta the tolerance of the derivative value
 * @param MaxIter the maximum number of iterations
 * 
 * @code {.cc}
 * NewtonSolver<SolverFunction> newton(f, x0, eps, delta, MaxIter);
 * @endcode
 */
template <class F>
NewtonSolver<F>::NewtonSolver(F f, double x0, double eps, double delta, int MaxIter): 
    f(f), x0(x0), eps(eps), delta(delta), MaxIter(MaxIter), iter(0), x(x0)
{}

/**
 * @brief Solve the equation.
 * 
 * @tparam F the function class
 * @return double the root
 * 
 * @code {.cc}
 * double root = newton.solve();
 * @endcode
 */
template <class F>
double NewtonSolver<F>::solve()
{
    double fx = f(x);
    double dfx = f.derivative(x);
    while (std::abs(fx) > eps && std::abs(dfx) > delta && iter < MaxIter)
    {
        x = x - fx / dfx;
        fx = f(x);
        dfx = f.derivative(x);
        iter++;
    }
    return x;
}

template class NewtonSolver<SolverFunction>;
template class NewtonSolver<Polynomial>;