#include "Polynomial.hpp"
#include "Lagrange.hpp"

Lagrange::Lagrange(const std::vector<double> &xData, const std::vector<double> &yData)
{
    x_lists = xData;
    y_lists = yData;
    m_poly = interpolate(xData, yData);
}

Polynomial Lagrange::interpolate(const std::vector<double> &xData, const std::vector<double> &yData)
{
    int n = xData.size();
    Polynomial result;
    for (auto i = 0; i < n; i++)
    {
        Polynomial term({yData[i]});
        for (auto j = 0; j < n; j++)
        {
            if (j == i)
            {
                continue;
            }
            Polynomial temp({-xData[j], 1});
            term = term * temp / (xData[i] - xData[j]);
        }
        result = result + term;
    }
    return result;
}

double Lagrange::operator()(double x) const
{
    return m_poly(x);
}

int Lagrange::degree() const
{
    return m_poly.get_degree();
}

std::ostream &operator<<(std::ostream &os, const Lagrange &lagrange)
{
    os << lagrange.m_poly;
    return os;
}