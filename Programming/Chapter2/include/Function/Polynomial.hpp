#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <iostream>
#include <vector>
#include <cmath>

class Polynomial
{
public:
    Polynomial();
    Polynomial(const std::vector<double> &coeffs);
    Polynomial(const Polynomial &rhs);
    Polynomial &operator=(const Polynomial &rhs);
    double operator()(double x) const;
    int get_degree() const;
    Polynomial operator+(const Polynomial &rhs) const;
    Polynomial operator-(const Polynomial &rhs) const;
    Polynomial operator*(const Polynomial &rhs) const;
    Polynomial operator*(double scalar) const;
    Polynomial operator/(double scalar) const;
    Polynomial derivative() const;
    Polynomial integral() const;
    double integrate(double a, double b) const;
    friend std::ostream &operator<<(std::ostream &os, const Polynomial &poly);

    std::vector<double> Get_all_roots() const;
private:
    std::vector<double> coefficients;
    int degree;
};

#endif