/**
 * @file SolverFunction.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Find the root of the function.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */


#ifndef __SOLVERFUNCTION_HPP__
#define __SOLVERFUNCTION_HPP__

#include "Function.hpp"
class SolverFunction : public Function
{
public:
    SolverFunction(double (*f)(double), double (*df)(double)) : f(f), df(df) {}
    SolverFunction(double (*f)(double)) : f(f), df(nullptr) {}

    inline double operator()(double x) const override { return f(x);} 
    inline double derivative(double x) const override
    { 
        if (df != nullptr) 
            return df(x);
        else
        {
            /**<< use \f$f^{\prime} (x_0) = \frac{f(x_0 + \Delta x) - f(x_0 - \Delta x)}{2\Delta x}\f$ */
            double delta_x = 1e-4;
            return (f(x+delta_x) + f(x-delta_x) - 2.0*f(x)) / (2.0 * delta_x);
        }
    }

    inline double integral(double a, double b) const override
    {
        double delta_x = 1e-3;
        double result = 0.0;
        for (double x = a; x < b; x += delta_x)
        {
            result += f(x) * delta_x;
        }
        return 0;
    }

private:
    double (*f)(double); // function pointer
    double (*df)(double); // derivative function pointer
};

#endif