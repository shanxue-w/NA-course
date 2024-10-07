#ifndef __POLYINTERPOLATION_HPP__
#define __POLYINTERPOLATION_HPP__

#include <iostream>
#include <vector>
#include "Polynomial.hpp"
// 虚函数，这是父类
class PolyInterpolation
{
public:
    virtual Polynomial interpolate(const std::vector<double> &xData, const std::vector<double> &yData) = 0;
    virtual double operator()(double x) const = 0;
    virtual int degree() const = 0;
    friend std::ostream &operator<<(std::ostream &os, const PolyInterpolation &polyInterpolation);
    virtual ~PolyInterpolation() {}
private:
    Polynomial m_poly;
};

#endif