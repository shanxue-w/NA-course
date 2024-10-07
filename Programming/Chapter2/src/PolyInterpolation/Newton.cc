#include "Newton.hpp"
#include "Polynomial.hpp"

template <class Poly>
Newton<Poly>::Newton()
{
    x_lists = {};
    y_lists = {};
    m_poly = Poly();
    divided_diff = {};
}

template <class Poly>
Newton<Poly>::Newton(const std::vector<double> &xData, const std::vector<double> &yData)
{
    x_lists = xData;
    y_lists = yData;
    m_poly = interpolate(xData, yData);
}

template <class Poly>
Poly Newton<Poly>::interpolate(const std::vector<double> &xData, const std::vector<double> &yData)
{
    int n = xData.size();
    divided_diff.resize(n);
    for (int i = 0; i < n; i++)
    {
        divided_diff[i].resize(i + 1);
        divided_diff[i][0] = yData[i];
    }
    for (int j = 1; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            divided_diff[i][j] = (divided_diff[i][j - 1] - divided_diff[i - 1][j - 1]) / (xData[i] - xData[i - j]);
        }
    }
    Poly result;
    for (int i = 0; i < n; i++)
    {
        Poly tmp({divided_diff[i][i]});
        for (int j = 0; j < i; j++)
        {
            tmp = tmp * Poly({-xData[j], 1});
        }
        result = result + tmp;
    }
    return result;
}

template <class Poly>
Poly Newton<Poly>::add_point(double x, double y)
{
    x_lists.push_back(x);
    y_lists.push_back(y);
    int n = x_lists.size();
    divided_diff.push_back({});
    divided_diff[n - 1].resize(n);
    for (int j = 1; j < n; j++)
    {
        divided_diff[n - 1][j] = (divided_diff[n - 1][j - 1] - divided_diff[n - 2][j - 1]) / (x - x_lists[n - 1 - j]);
    }
    Poly tmp({divided_diff[n - 1][n - 1]});
    for (int j = 0; j < n - 1; j++)
    {
        tmp = tmp * Poly({-x_lists[j], 1});
    }
    return m_poly + tmp;
}

template <class Poly>
double Newton<Poly>::operator()(double x) const
{
    return m_poly(x);
}

template <class Poly>
int Newton<Poly>::degree() const
{
    return m_poly.get_degree();
}

// template <class P>
// std::ostream &operator<<(std::ostream &os, const Newton<P> &newton)
// {
//     os << newton.m_poly;
//     return os;
// }

template <class Poly>
double Newton<Poly>::derivative(double x, int n) const
{
    Poly poly = m_poly;
    for (int i = 0; i < n; i++)
    {
        poly = poly.derivative();
    }
    return poly(x);
}

template class Newton<Polynomial>;