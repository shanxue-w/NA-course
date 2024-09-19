/**
 * @file EquationSolver.h
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2024-09-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef EquationSolver_H
#define EquationSolver_H

#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "Function.hpp"

class EquationSolver
{
public:
    virtual double solve() = 0;

    virtual ~EquationSolver() {}
};

#endif