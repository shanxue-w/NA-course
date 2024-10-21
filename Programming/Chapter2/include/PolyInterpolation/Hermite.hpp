#ifndef __HERMITE_HPP__
#define __HERMITE_HPP__


#include "PolyInterpolation.hpp"
#include "Polynomial.hpp"
#include "NewtonPoly.hpp"
#include <tuple>


class Hermite : public PolyInterpolation
{
public:
    Hermite();
    // 三元组 (x,y,n) n 表示 y是f的第n阶导数
    Hermite(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<int> &nData);
    NewtonPoly interpolate(std::vector<std::tuple<double, double, int>> &data_lists);
    NewtonPoly add_point(double x, double y, int n);
    double operator()(double x) const override;
    NewtonPoly get_polynomial(int n=0) const;
    double derivative(double x) const;
    double integral(double a, double b) const;
    int degree() const override;
    NewtonPoly get_polynomial();
    Polynomial Convert_to_Polynomial() const;

    friend std::ostream &operator<<(std::ostream &os, const Hermite &hermite)
    {
        os << hermite.m_poly;
        return os;
    }
private:
    std::vector<std::tuple<double, double, int>> data; // 三元组 (x,y,n) n 表示 y是f的第n阶导数
    std::vector<std::vector<double>> divided_diff; // divided difference table
    std::vector<double> x_lists;
    std::vector<double> y_lists;
    NewtonPoly m_poly;
};

#endif