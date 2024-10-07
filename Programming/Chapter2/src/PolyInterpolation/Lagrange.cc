#include "Lagrange.hpp"
#include "Polynomial.hpp"

template <class Poly>
Lagrange<Poly>::Lagrange()
{
    x_lists = {};
    y_lists = {};
    m_poly = Poly();
}

template <class Poly>
Lagrange<Poly>::Lagrange(const std::vector<double> &xData, const std::vector<double> &yData)
{
    x_lists = xData;
    y_lists = yData;
    m_poly = interpolate(xData, yData);
}

template <class Poly>
Poly Lagrange<Poly>::interpolate(const std::vector<double> &xData, const std::vector<double> &yData)
{
    int n = xData.size();
    Poly result;
    for (auto i = 0; i < n; i++)
    {
        Poly term({yData[i]});
        for (auto j = 0; j < n; j++)
        {
            if (j == i)
            {
                continue;
            }
            Poly temp({-xData[j], 1});
            term = term * temp / (xData[i] - xData[j]);
        }
        result = result + term;
    }
    return result;
}

template <class Poly>
double Lagrange<Poly>::operator()(double x) const
{
    return m_poly(x);
}

template <class Poly>
int Lagrange<Poly>::degree() const
{
    return m_poly.get_degree();
}

// template <class P>
// std::ostream &operator<<(std::ostream &os, const Lagrange<P> &lagrange)
// {
//     os << lagrange.m_poly;
//     return os;
// }

template <class Poly>
double Lagrange<Poly>::derivative(double x, int n) const
{
    Poly poly = m_poly;
    for (int i = 0; i < n; i++)
    {
        poly = poly.derivative();
    }
    return poly(x);
}

template <class Poly>
double Lagrange<Poly>::integral(double a, double b) const
{
    Poly poly = m_poly;
    poly = poly.integral();
    return poly(b) - poly(a);
}

template class Lagrange<Polynomial>;