#ifndef __NEWTON_HPP__
#define __NEWTON_HPP__

#include "PolyInterpolation.hpp"
#include <iostream>
class Newton : public PolyInterpolation
{
public:
    Newton() {}
    Newton(const std::vector<double> &xData, const std::vector<double> &yData);
    Polynomial interpolate(const std::vector<double> &xData, const std::vector<double> &yData) override;
    Polynomial add_point(double x, double y);
    double operator()(double x) const override;
    int degree() const override;
    friend std::ostream &operator<<(std::ostream &os, const Newton &newton);
    ~Newton() {}
};

#endif