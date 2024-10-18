#include "NewtonSolver.hpp"
#include "SolverFunction.hpp"
#include "Polynomial.hpp"

template <class F>
NewtonSolver<F>::NewtonSolver(F f, double x0, double eps, double delta, int MaxIter): 
    f(f), x0(x0), eps(eps), delta(delta), MaxIter(MaxIter), iter(0), x(x0)
{}

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