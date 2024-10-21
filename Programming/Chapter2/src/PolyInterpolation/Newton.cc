#include "Newton.hpp"
#include "Polynomial.hpp"
#include "NewtonPoly.hpp"

Newton::Newton()
{
    x_lists = {};
    y_lists = {};
    m_poly = NewtonPoly();
    divided_diff = {};
}


Newton::Newton(const std::vector<double> &xData, const std::vector<double> &yData)
{
    x_lists = xData;
    y_lists = yData;
    m_poly = interpolate(xData, yData);
}

NewtonPoly Newton::interpolate(const std::vector<double> &xData, const std::vector<double> &yData)
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
    std::vector<double> w_lists;
    for (int i = 0; i < n; i++)
    {
        w_lists.push_back(divided_diff[i][i]);
    }
    NewtonPoly newtonPoly(xData, yData, w_lists);
    return newtonPoly;
}


NewtonPoly Newton::add_point(double x, double y)
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
    std::vector<double> w_lists;
    for (int i = 0; i < n; i++)
    {
        w_lists.push_back(divided_diff[i][i]);
    }
    NewtonPoly newtonPoly(x_lists, y_lists, w_lists);
    return newtonPoly;
}


double Newton::operator()(double x) const
{
    return m_poly(x);
}


int Newton::degree() const
{
    return m_poly.get_degree();
}

// template <class P>
// std::ostream &operator<<(std::ostream &os, const Newton<P> &newton)
// {
//     os << newton.m_poly;
//     return os;
// }


double Newton::derivative(double x) const
{
    return m_poly.derivative(x);
}


double Newton::integral(double a, double b) const
{
    return m_poly.integral(a, b);
}

NewtonPoly Newton::getNewtonPoly() const
{
    return m_poly;
}

Polynomial Newton::Convert_to_Polynomial() const
{
    return m_poly.Convert_to_Polynomial();
}