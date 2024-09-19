/**
 * @file BisectionSolver.h
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BisectionSolver_H
#define BisectionSolver_H

#include "EquationSolver.h"

class BisectionSolver : public EquationSolver
{
public:
    BisectionSolver (double (*f)(double), double a, double b); // 构造函数，传入函数指针和区间
    BisectionSolver (double (*f)(double), double a, double b, int MaxIter); // 构造函数，传入函数指针和区间
    double solve() override; // 重载求解函数
    double getRoot() const { return x; } // 获取解
    int    getIter() const { return iter; } // 获取迭代次数
private:
    double (*f)(double); // 函数指针
    double a, b; // 区间
    double eps;
    int MaxIter;
    int    iter; // 迭代次数
    double x; // 解
};


#endif 