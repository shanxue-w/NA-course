/**
 * @file Newton.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The Newton class, providing the basic operations of Newton Interpolation.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __NEWTON_HPP__
#define __NEWTON_HPP__

#include "PolyInterpolation.hpp"
#include "NewtonPoly.hpp"
#include <iostream>



class Newton : public PolyInterpolation
{
public:
    Newton();
    Newton(const std::vector<double> &xData, const std::vector<double> &yData);
    Newton(const std::vector<double> &xData, FuncPtr f);
    NewtonPoly interpolate(const std::vector<double> &xData, const std::vector<double> &yData);
    NewtonPoly add_point(double x, double y);
    double operator()(double x) const override;
    double derivative(double x) const;
    double integral(double a, double b) const;
    int degree() const override;
    NewtonPoly getNewtonPoly() const;
    Polynomial Convert_to_Polynomial() const;
    
    friend std::ostream &operator<<(std::ostream &os, const Newton &newton)
    {
        os << newton.m_poly;
        return os;
    }
    ~Newton() {}

private:
    std::vector<double> x_lists;
    std::vector<double> y_lists;
    // 函数指针,double
    FuncPtr f;
    std::vector<std::vector<double>> divided_diff; // divided difference table
    NewtonPoly m_poly;
};

#endif