/**
 * @file LagPoly.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The class for Lagrange Polynomial, which is not used in the homework.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __LAGPOLY_HPP__
#define __LAGPOLY_HPP__

#include "Function.hpp"
#include <vector>
#include <iostream>

class LagPoly : public Function
{
public:
    LagPoly();
    LagPoly(const std::vector<double> xData, const std::vector<double> yData);
    LagPoly(const std::vector<double> xData, const std::vector<double> yData, const std::vector<double> wData);
    LagPoly &operator=(const LagPoly &rhs);
    double operator()(double x) const override;
    double derivative(double x) const override;
    double integral(double a, double b) const override;
    int get_degree() const;
    friend std::ostream &operator<<(std::ostream &os, const LagPoly &lagpoly);

private:
    std::vector<double> x_lists;
    std::vector<double> y_lists;
    std::vector<double> w_lists;
};

#endif