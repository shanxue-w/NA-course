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