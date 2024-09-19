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

#include "EquationSolver.h"

class NewtonSolver : public EquationSolver
{
public:
    NewtonSolver (double (*f)(double), double (*df)(double), double x0); // 构造函数，传入函数指针和初始值
    NewtonSolver (double (*f)(double), double (*df)(double), double x0, int MaxIter); // 构造函数，传入函数指针和初始值
    double solve() override; // 重载求解函数
    double getRoot() const { return x; } // 获取解
    int    getIter() const { return iter; } // 获取迭代次数
private:
    double (*f)(double); // 函数指针
    double (*df)(double); // 函数导数指针
    double x0; // 初始值
    double eps;
    int MaxIter;
    int    iter; // 迭代次数
    double x; // 解
};

#endif