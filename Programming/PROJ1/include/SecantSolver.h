/**
 * @file SecantSolver.h
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SecantSolver_H
#define SecantSolver_H

#include "EquationSolver.h"

class SecantSolver : public EquationSolver
{
public:
    SecantSolver (double (*f)(double), double x0, double x1); // 构造函数，传入函数指针和初始值
    SecantSolver (double (*f)(double), double x0, double x1, int MaxIter); // 构造函数，传入函数指针和初始值
    double solve() override; // 重载求解函数
    double getRoot() const { return x; } // 获取解
    int    getIter() const { return iter; } // 获取迭代次数
private:
    double (*f)(double); // 函数指针
    double x0, x1; // 初始值
    double eps;
    int MaxIter;
    int    iter; // 迭代次数
    double x; // 解
};

#endif