#ifndef __LANGRANGE_HPP__
#define __LANGRANGE_HPP__

#include "PolyInterpolation.hpp"
#include <iostream>
class Lagrange : public PolyInterpolation 
{
public:
    Lagrange() {}
    Lagrange(const std::vector<double> &xData, const std::vector<double> &yData);
    Polynomial interpolate(const std::vector<double> &xData, const std::vector<double> &yData) override;
    double operator()(double x) const override;
    int degree() const override;
    friend std::ostream &operator<<(std::ostream &os, const Lagrange &lagrange);
    ~Lagrange() {}

private:
    std::vector<double> x_lists;
    std::vector<double> y_lists;
    Polynomial          m_poly;
};

#endif