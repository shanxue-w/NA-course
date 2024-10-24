/**
 * @file NewtonSolver.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The NewtonSolver class, providing the basic operations of Newton Solver.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __NEWTONSOLVER_HPP__
#define __NEWTONSOLVER_HPP__

#include "EquationSolver.hpp"
#include "SolverFunction.hpp"

template <class F>
class NewtonSolver : public EquationSolver
{
public:
    NewtonSolver (F f, double x0, double eps=1e-12, double delta=1e-12, int MaxIter=100); // 构造函数
    double solve() override; // 重载求解函数
    inline double getRoot() const { return x; } // 获取解
    inline int    getIter() const { return iter; } // 获取迭代次数
    inline double operator()(double x) const { return f(x); } // 重载运算符
private:
    F f; // 函数指针
    double x0; // 初始值
    double eps;
    double delta;
    int MaxIter;
    int    iter; // 迭代次数
    double x; // 解
};

#endif