/**
 * @file Function.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Define the function class
 * @version 0.1
 * @date 2024-09-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FUNCTION_HPP
#define FUNCTION_HPP

class Function
{
public:
    Function(double (*f)(double), double (*df)(double)) : f(f), df(df) {}
    Function(double (*f)(double)) : f(f), df(nullptr) {}

    double operator()(double x) { return f(x); }
    double derivative(double x) { return df(x); }

private:
    double (*f)(double); // function pointer
    double (*df)(double); // derivative function pointer
};

#endif