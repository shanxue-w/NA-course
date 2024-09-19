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

#include "EquationSolver.hpp"

class BisectionSolver : public EquationSolver
{
public:
    BisectionSolver (Function f, double a, double b, double eps=1e-12, double delta=1e-9, int MaxIter=100); // 构造函数
    double solve() override; // 重载求解函数
    double getRoot() const { return x; } // 获取解
    int    getIter() const { return iter; } // 获取迭代次数
private:
    Function f; /*<< The function */
    double   a, b; /*<< */
    double   eps;
    double   delta;
    int      MaxIter;
    int      iter; // 迭代次数
    double   x; // 解
};


#endif 