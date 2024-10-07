#ifndef __NEWTON_HPP__
#define __NEWTON_HPP__

#include "PolyInterpolation.hpp"
#include <iostream>

template <class Poly>
class Newton : public PolyInterpolation<Poly>
{
public:
    Newton();
    Newton(const std::vector<double> &xData, const std::vector<double> &yData);
    Poly interpolate(const std::vector<double> &xData, const std::vector<double> &yData) override;
    Poly add_point(double x, double y);
    double operator()(double x) const override;
    double derivative(double x, int n) const;
    double integral(double a, double b) const;
    int degree() const override;
    
    template <class P>
    friend std::ostream &operator<<(std::ostream &os, const Newton<P> &newton)
    {
        os << newton.m_poly;
        return os;
    }
    ~Newton() {}

private:
    std::vector<double> x_lists;
    std::vector<double> y_lists;
    Poly          m_poly;
    std::vector<std::vector<double>> divided_diff; // divided difference table
};

#endif