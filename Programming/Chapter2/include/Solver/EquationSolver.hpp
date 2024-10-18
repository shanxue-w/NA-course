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