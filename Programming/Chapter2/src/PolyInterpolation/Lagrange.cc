#include "LagPoly.hpp"
#include "Lagrange.hpp"

Lagrange::Lagrange() : poly() {}

Lagrange::Lagrange(const std::vector<double> &xData, const std::vector<double> &yData) : poly(xData, yData) {}

Lagrange::Lagrange(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &wData) : 
poly(xData, yData, wData) 
{}

double Lagrange::operator()(double x) const
{
    return poly(x);
}

int Lagrange::degree() const
{
    return poly.get_degree();
}

LagPoly Lagrange::get_polynomial() const
{
    return poly;
}

std::ostream &operator<<(std::ostream &os, const Lagrange &lagrange)
{
    os << lagrange.poly;
    return os;
}
