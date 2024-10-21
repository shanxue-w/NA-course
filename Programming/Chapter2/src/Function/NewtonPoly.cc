#include "NewtonPoly.hpp"
#include "NewtonSolver.hpp"
#include <numeric>

NewtonPoly::NewtonPoly() : x_lists({0}), y_lists({0}), w_lists({0}) {}

NewtonPoly::NewtonPoly(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &wData) : 
x_lists(xData), y_lists(yData), w_lists(wData) 
{}

NewtonPoly::NewtonPoly(const NewtonPoly &rhs) : x_lists(rhs.x_lists), y_lists(rhs.y_lists), w_lists(rhs.w_lists) {}

NewtonPoly &NewtonPoly::operator=(const NewtonPoly &rhs)
{
    if (this != &rhs)
    {
        x_lists = rhs.x_lists;
        y_lists = rhs.y_lists;
        w_lists = rhs.w_lists;
    }
    return *this;
}

double NewtonPoly::operator()(double x) const
{
    double result = 0.0;
    int n = x_lists.size();
    for (auto i= n-1; i >= 0; i--)
    {
        result = result * (x - x_lists[i]) + w_lists[i];
    }
    return result;
}

double NewtonPoly::derivative(double x) const
{
    double epison = 1e-4;
    return (this->operator()(x+epison) - this->operator()(x-epison)) / (2*epison);
}

double NewtonPoly::integral(double a, double b) const
{
    double delta_x = 1e-3;
    double result = 0.0;
    for (double x = a; x < b; x += delta_x)
    {
        result += this->operator()(x) * delta_x;
    }
    return result;
}

int NewtonPoly::get_degree() const
{
    return x_lists.size()-1;
}

std::ostream &operator<<(std::ostream &os, const NewtonPoly &newtonpoly)
{
    int n = newtonpoly.x_lists.size();
    for (auto i=0; i<n; i++)
    {
        if (i == 0)
        {
            os << newtonpoly.w_lists[i];
        }
        else
        {
            if (newtonpoly.w_lists[i] > 0)
            {
                os << " + " << newtonpoly.w_lists[i];
            }
            else
            {
                os << " - " << -newtonpoly.w_lists[i];
            }
        }
        for (auto j=0; j<i; j++)
        {
            if (newtonpoly.x_lists[j] < 0)
            {
                os << "(x + " << -newtonpoly.x_lists[j] << ")";
            }
            else
            {
                os << "(x - " << newtonpoly.x_lists[j] << ")";
            }
        }
    }
    return os;
}

Polynomial NewtonPoly::Convert_to_Polynomial() const
{
    int n = x_lists.size();
    Polynomial poly({w_lists[0]});
    for (auto i=1; i<n; i++)
    {
        Polynomial temp({w_lists[i]});
        for (auto j=0; j<i; j++)
        {
            Polynomial factor({-x_lists[j], 1});
            temp = temp * factor;
        }
        poly = poly + temp;
    }
    return poly;
}