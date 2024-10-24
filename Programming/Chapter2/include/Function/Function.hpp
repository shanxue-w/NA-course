/**
 * @file Function.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The abstract class of function.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __FUNCTION_HPP__
#define __FUNCTION_HPP__
using FuncPtr = double (*)(double);
class Function
{
public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const = 0;
    virtual double integral(double a, double b) const = 0;
    virtual ~Function() {}
private:
    
};

#endif