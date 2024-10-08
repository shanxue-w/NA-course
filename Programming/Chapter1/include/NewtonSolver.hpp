/**
 * @file NewtonSolver.h
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NewtonSolver_H
#define NewtonSolver_H

#include "EquationSolver.hpp"

class NewtonSolver : public EquationSolver
{
public:
    NewtonSolver (Function f, double x0, double eps=1e-12, double delta=1e-9, int MaxIter=100); // 构造函数
    double solve() override; // 重载求解函数
    inline double getRoot() const { return x; } // 获取解
    inline int    getIter() const { return iter; } // 获取迭代次数
private:
    Function f; // 函数指针
    double x0; // 初始值
    double eps;
    double delta;
    int MaxIter;
    int    iter; // 迭代次数
    double x; // 解
};

#endif