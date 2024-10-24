/**
 * @file Lagrange.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The Lagrange class, providing the basic operations of Lagrange Interpolation.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __LANGRANGE_HPP__
#define __LANGRANGE_HPP__

#include "PolyInterpolation.hpp"
#include "LagPoly.hpp"
#include <iostream>

class Lagrange : public PolyInterpolation
{
public:
    Lagrange();
    Lagrange(const std::vector<double> &xData, const std::vector<double> &yData);
    Lagrange(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &wData);
    double operator()(double x) const override;
    int degree() const override;
    LagPoly get_polynomial() const;
    friend std::ostream &operator<<(std::ostream &os, const Lagrange &lagrange);
private:
    LagPoly poly;
};

#endif