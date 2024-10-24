/**
 * @file EquationSolver.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The EquationSolver class, providing the basic operations of EquationSolver.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __EQUATIONSOLVER_HPP__
#define __EQUATIONSOLVER_HPP__

#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "SolverFunction.hpp"

class EquationSolver
{
public:
    virtual double solve() = 0;

    virtual ~EquationSolver() {}
};

#endif