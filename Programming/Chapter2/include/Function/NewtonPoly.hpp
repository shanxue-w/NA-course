/**
 * @file NewtonPoly.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The class for Newton Polynomial, providing the basic operations of Newton Polynomial.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

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