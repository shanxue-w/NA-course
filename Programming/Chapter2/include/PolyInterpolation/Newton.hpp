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
    std::vector<std::vector<double>> divided_diff; // divided difference table
    NewtonPoly m_poly;
};

#endif