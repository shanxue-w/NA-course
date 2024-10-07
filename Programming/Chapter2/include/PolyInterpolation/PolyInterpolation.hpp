#ifndef __POLYINTERPOLATION_HPP__
#define __POLYINTERPOLATION_HPP__

#include "Polynomial.hpp"
#include <iostream>
#include <vector>
// 虚函数，这是父类

template <class Poly>
class PolyInterpolation
{
public:
    virtual Poly interpolate(const std::vector<double> &xData, const std::vector<double> &yData) = 0;
    virtual double operator()(double x) const = 0;
    virtual int degree() const = 0;
    // friend std::ostream &operator<< <> (std::ostream &os, const PolyInterpolation<Poly> &polyInterpolation);
    virtual ~PolyInterpolation() {}
private:
    Poly m_poly;
};

#endif