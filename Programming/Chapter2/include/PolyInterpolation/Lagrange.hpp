#ifndef __LANGRANGE_HPP__
#define __LANGRANGE_HPP__

#include "PolyInterpolation.hpp"
#include <iostream>

template <class Poly>
class Lagrange : public PolyInterpolation<Poly>
{
public:
    Lagrange();
    Lagrange(const std::vector<double> &xData, const std::vector<double> &yData);
    Poly interpolate(const std::vector<double> &xData, const std::vector<double> &yData) override;
    double operator()(double x) const override;
    double derivative(double x, int n) const;
    double integral(double a, double b) const;
    int degree() const override;

    template <class P>
    friend std::ostream &operator<<(std::ostream &os, const Lagrange<P> &lagrange)
    {
        os << lagrange.m_poly;
        return os;
    }
    ~Lagrange() {}

private:
    std::vector<double> x_lists;
    std::vector<double> y_lists;
    Poly          m_poly;
};

#endif