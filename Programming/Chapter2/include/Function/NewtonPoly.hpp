#ifndef __NEWTONPOLY_HPP__
#define __NEWTONPOLY_HPP__

#include "Function.hpp"
#include "Polynomial.hpp"
#include <vector>
#include <iostream>


class NewtonPoly : public Function
{
public:
    NewtonPoly();
    NewtonPoly(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &wData);
    NewtonPoly(const NewtonPoly &rhs);
    NewtonPoly &operator=(const NewtonPoly &rhs);
    double operator()(double x) const override;
    double derivative(double x) const override;
    double integral(double a, double b) const override;
    int get_degree() const;
    friend std::ostream &operator<<(std::ostream &os, const NewtonPoly &newtonpoly);
    Polynomial Convert_to_Polynomial() const;
    ~NewtonPoly() override = default;

private:
    std::vector<double> x_lists;
    std::vector<double> y_lists;
    std::vector<double> w_lists;
};

#endif