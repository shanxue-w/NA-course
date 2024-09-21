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

#include "EquationSolver.hpp"

class SecantSolver : public EquationSolver
{
public:
    SecantSolver (Function f, double x0, double x1, double eps=1e-12, double delta=1e-9, int MaxIter=100); 
    double solve() override; 
    inline double getRoot() const { return x; } 
    inline int    getIter() const { return iter; }
private:
    Function f;
    double x0, x1; // 初始值
    double eps;
    double delta;
    int MaxIter;
    int    iter; // 迭代次数
    double x; // 解
};

#endif