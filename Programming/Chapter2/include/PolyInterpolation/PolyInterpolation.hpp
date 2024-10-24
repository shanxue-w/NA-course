/**
 * @file PolyInterpolation.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The PolyInterpolation class, providing the basic operations of PolyInterpolation.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __POLYINTERPOLATION_HPP__
#define __POLYINTERPOLATION_HPP__

#include "Polynomial.hpp"
// #include "LagPoly.hpp"
#include <iostream>
#include <vector>
// 虚函数，这是父类

class PolyInterpolation
{
public:
    // PolyInterpolation(const std::vector<double> &xData, const std::vector<double> &yData);
    virtual double operator()(double x) const = 0;
    virtual int degree() const = 0;
    // friend std::ostream &operator<< <> (std::ostream &os, const PolyInterpolation<Poly> &polyInterpolation);
    virtual ~PolyInterpolation() {}
};

#endif