#ifndef __HERMITE_HPP__
#define __HERMITE_HPP__


#include "PolyInterpolation.hpp"
#include <tuple>

template <class Poly>
class Hermite : public PolyInterpolation<Poly>
{
public:
    Hermite();
    // 三元组 (x,y,n) n 表示 y是f的第n阶导数
    Hermite(const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<int> &nData);
    Poly interpolate(const std::vector<double> &xData, const std::vector<double> &yData) override;
    Poly add_point(double x, double y, int n);
    double operator()(double x) const override;
    Poly get_polynomial(int n=0) const;
    double derivative(double x, int n) const;
    double integral(double a, double b) const;
    int degree() const override;

    template <class P>
    friend std::ostream &operator<<(std::ostream &os, const Hermite<P> &hermite)
    {
        os << hermite.m_poly;
        return os;
    }
private:
    std::vector<std::tuple<double, double, int>> data; // 三元组 (x,y,n) n 表示 y是f的第n阶导数
    Poly          m_poly;
    std::vector<std::vector<double>> divided_diff; // divided difference table
    std::vector<double> x_lists;
    std::vector<double> y_lists;
};

#endif