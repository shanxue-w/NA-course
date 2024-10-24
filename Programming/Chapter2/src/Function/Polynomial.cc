/**
 * @file Polynomial.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief This is the implementation of the Polynomial class.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "Polynomial.hpp"
#include "NewtonSolver.hpp"

/**
 * @brief Construct a new Polynomial:: Polynomial object
 * 
 * @code {.cc}
 * Polynomial poly;
 * @endcode
 * 
 */
Polynomial::Polynomial() : coefficients({0}), degree(0) {}

/**
 * @brief Construct a new Polynomial:: Polynomial object
 * 
 * @param coeffs the coefficients of the polynomial, base is
 * \f[
 * 1, x, x^2, x^3, \cdots, x^n
 * \f]
 * 
 * @code {.cc}
 * std::vector<double> coeffs = {1, 2, 3}; // 1 + 2x + 3x^2
 * Polynomial poly(coeffs);
 * @endcode
 * 
 */
Polynomial::Polynomial(const std::vector<double> &coeffs) : coefficients(coeffs), degree(coeffs.size() - 1) {}

/**
 * @brief Construct a new Polynomial:: Polynomial object
 * 
 * @param rhs The right hand side of the assignment operator
 * 
 * @code {.cc}
 * Polynomial poly1;
 * Polynomial poly2(poly1);
 * @endcode
 */
Polynomial::Polynomial(const Polynomial &rhs) : coefficients(rhs.coefficients), degree(rhs.degree) {}

/**
 * @brief Construct a new Polynomial:: Polynomial object
 * 
 * @param rhs The right hand side of the assignment operator
 * @return Polynomial& 
 * 
 * @code {.cc}
 * Polynomial poly1;
 * Polynomial poly2;
 * poly1 = poly2;
 * @endcode
 */
Polynomial &Polynomial::operator=(const Polynomial &rhs)
{
    if (this != &rhs)
    {
        coefficients = rhs.coefficients;
        degree = rhs.degree;
    }
    return *this;
}

/**
 * @brief Use the Horner's method to calculate the value of the Polynomial.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * double y = poly(1.0);
 * @endcode
 */
double Polynomial::operator()(double x) const
{
    double result = 0;
    for (int i = degree; i >= 0; i--)
    {
        result = result * x + coefficients[i];
    }
    return result;
}

/**
 * @brief Calculate the derivative of the Polynomial.
 * 
 * @param x the value of x
 * @return double 
 * 
 * @code {.cc}
 * double y = poly.derivative(1.0);
 * @endcode
 */
double Polynomial::derivative(double x) const
{
    double result = 0;
    for (int i = degree; i >= 1; i--)
    {
        result = result * x + coefficients[i] * i;
    }
    return result;
}

/**
 * @brief Calculate the integral of the Polynomial.
 * 
 * @param a the lower bound
 * @param b the upper bound
 * @return double 
 * 
 * @code {.cc}
 * double y = poly.integral(0.0, 1.0);
 * @endcode
 */
double Polynomial::integral(double a, double b) const
{
    Polynomial integral = this->integral();
    return integral(b) - integral(a);
}

/**
 * @brief Get the degree of the Polynomial.
 * 
 * @return int 
 * 
 * @code {.cc}
 * int n = poly.get_degree();
 * @endcode
 */
int Polynomial::get_degree() const
{
    return degree;
}

/**
 * @brief overload the operator+ to add two Polynomials.
 * 
 * @param rhs the right hand side of the operator+
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2({3, 2, 1}); // 3 + 2x + x^2
 * Polynomial poly3 = poly1 + poly2; // 4 + 4x + 4x^2
 * @endcode
 */
Polynomial Polynomial::operator+(const Polynomial &rhs) const
{
    std::vector<double> result(std::max(degree, rhs.degree) + 1, 0);
    for (int i = 0; i <= degree; i++)
    {
        result[i] += coefficients[i];
    }
    for (int i = 0; i <= rhs.degree; i++)
    {
        result[i] += rhs.coefficients[i];
    }
    return Polynomial(result);
}

/**
 * @brief overload the operator- to subtract two Polynomials.
 * 
 * @param rhs the right hand side of the operator-
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2({3, 2, 1}); // 3 + 2x + x^2
 * Polynomial poly3 = poly1 - poly2; // -2 + 0x + 2x^2
 * @endcode
 */
Polynomial Polynomial::operator-(const Polynomial &rhs) const
{
    std::vector<double> result(std::max(degree, rhs.degree) + 1, 0);
    for (int i = 0; i <= degree; i++)
    {
        result[i] += coefficients[i];
    }
    for (int i = 0; i <= rhs.degree; i++)
    {
        result[i] -= rhs.coefficients[i];
    }
    return Polynomial(result);
}

/**
 * @brief overload the operator* to multiply two Polynomials.
 * 
 * @param rhs the right hand side of the operator*
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2({3, 2, 1}); // 3 + 2x + x^2
 * Polynomial poly3 = poly1 * poly2; // 3 + 8x + 14x^2 + 8x^3 + 3x^4
 * @endcode
 */
Polynomial Polynomial::operator*(const Polynomial &rhs) const
{
    std::vector<double> result(degree + rhs.degree + 1, 0);
    for (int i = 0; i <= degree; i++)
    {
        for (int j = 0; j <= rhs.degree; j++)
        {
            result[i + j] += coefficients[i] * rhs.coefficients[j];
        }
    }
    return Polynomial(result);
}

/**
 * @brief overload the operator* to multiply the Polynomial by a scalar.
 * 
 * @param scalar the scalar
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2 = poly1 * 2; // 2 + 4x + 6x^2
 * @endcode
 */
Polynomial Polynomial::operator*(double scalar) const
{
    std::vector<double> result(coefficients);
    for (int i = 0; i <= degree; i++)
    {
        result[i] *= scalar;
    }
    return Polynomial(result);
}

/**
 * @brief overload the operator/ to divide the Polynomial by a scalar.
 * 
 * @param scalar the scalar
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2 = poly1 / 2; // 0.5 + x + 1.5x^2
 * @endcode
 */
Polynomial Polynomial::operator/(double scalar) const
{
    std::vector<double> result(coefficients);
    for (int i = 0; i <= degree; i++)
    {
        result[i] /= scalar;
    }
    return Polynomial(result);
}

/**
 * @brief Derivative of the Polynomial.
 * 
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2 = poly1.derivative(); // 2 + 6x
 * @endcode
 */
Polynomial Polynomial::derivative() const
{
    std::vector<double> result(std::max(degree - 1, 0) + 1, 0);
    for (int i = 1; i <= degree; i++)
    {
        result[i - 1] = coefficients[i] * i;
    }
    return Polynomial(result);
}

/**
 * @brief Integral of the Polynomial. The constant of integration is set to be 0.
 * 
 * @return Polynomial 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2 = poly1.integral(); // 0 + x + 1.5x^2
 * @endcode
 */
Polynomial Polynomial::integral() const
{
    std::vector<double> result(degree + 2, 0);
    for (int i = 0; i <= degree; i++)
    {
        result[i + 1] = coefficients[i] / (i + 1);
    }
    return Polynomial(result);
}

/**
 * @brief Calculate the integral of the Polynomial from a to b.
 * 
 * @param a the lower bound
 * @param b the upper bound
 * @return double 
 * 
 * @code {.cc}
 * double y = poly.integrate(0.0, 1.0);
 * @endcode
 */
double Polynomial::integrate(double a, double b) const
{
    Polynomial integral = this->integral();
    return integral(b) - integral(a);
}

/**
 * @brief Overload the operator<< to output the Polynomial.
 * 
 * @param os 
 * @param poly 
 * @return std::ostream& 
 * 
 * @code {.cc}
 * Polynomial poly({1, 2, 3}); // 1 + 2x + 3x^2
 * std::cout << poly << std::endl;
 * @endcode
 */
std::ostream &operator<<(std::ostream &os, const Polynomial &poly)
{
    int flag = 0;
    for (int i = poly.degree; i >= 0; i--)
    {
        if (std::abs(poly.coefficients[i]) < 1e-16)
        {
            continue;
        }
        if (i == poly.degree)
        {
            flag = 1;
            if (i == 1)
                os << poly.coefficients[i] << "x";
            else if (i == 0)
                os << poly.coefficients[i];
            else
                os << poly.coefficients[i] << "x^" << i;
        }
        else if (i == 1)
        {
            if (poly.coefficients[i] > 0)
            {
                if (flag == 0)
                {
                    flag = 1;
                    os << poly.coefficients[i] << "x";
                }
                else
                {
                    os << " + " << poly.coefficients[i] << "x";
                }
            }
            else if (std::abs(poly.coefficients[i]) < 1e-16)
            {
                continue;
            }
            else
            {
                if (flag == 0)
                {
                    flag = 1;
                    os << -poly.coefficients[i] << "x";
                }
                else
                {
                    os << " - " << -poly.coefficients[i] << "x";
                }
            }
        }
        else if (i == 0)
        {
            if (poly.coefficients[i] > 0)
            {
                if (flag == 0)
                {
                    flag = 1;
                    os << poly.coefficients[i];
                }
                else
                {
                    os << " + " << poly.coefficients[i];
                }
            }
            else if (std::abs(poly.coefficients[i]) < 1e-16)
            {
                continue;
            }
            else
            {
                if (flag == 0)
                {
                    flag = 1;
                    os << -poly.coefficients[i];
                }
                else
                {
                    os << " - " << -poly.coefficients[i];
                }
            }
        }
        else
        {
            if (poly.coefficients[i] > 0)
            {
                if (flag == 0)
                {
                    flag = 1;
                    os << poly.coefficients[i] << "x^" << i;
                }
                else
                {
                    os << " + " << poly.coefficients[i] << "x^" << i;
                }
            }
            else if (std::abs(poly.coefficients[i]) < 1e-16)
            {
                continue;
            }
            else
            {
                if (flag == 0)
                {
                    flag = 1;
                    os << -poly.coefficients[i] << "x^" << i;
                }
                else
                {
                    os << " - " << -poly.coefficients[i] << "x^" << i;
                }
            }
        }
    }
    return os;
}

/**
 * @brief Divide the Polynomial by another Polynomial.
 * 
 * @param rhs the right hand side of the operator/
 * @return std::pair<Polynomial, Polynomial> 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2({1, 1}); // 1 + x
 * std::pair<Polynomial, Polynomial> result = poly1.Divide(poly2);
 * Polynomial quotient = result.first;
 * Polynomial remainder = result.second;
 * // poly1 = poly2 * quotient + remainder
 * @endcode
 */
std::pair<Polynomial, Polynomial> Polynomial::Divide(const Polynomial &rhs) const
{
    std::vector<double> quotient(degree - rhs.degree + 1, 0);
    std::vector<double> remainder(coefficients);
    for (int i = degree - rhs.degree; i >= 0; i--)
    {
        quotient[i] = remainder[rhs.degree + i] / rhs.coefficients[rhs.degree];
        for (int j = 0; j <= rhs.degree; j++)
        {
            remainder[j + i] -= quotient[i] * rhs.coefficients[j];
        }
    }
    return std::make_pair(Polynomial(quotient), Polynomial(remainder));
}

/**
 * @brief overload the operator/ to divide the Polynomial by another Polynomial.
 * 
 * @param rhs the right hand side of the operator/
 * @return std::pair<Polynomial, Polynomial> 
 * 
 * @code {.cc}
 * Polynomial poly1({1, 2, 3}); // 1 + 2x + 3x^2
 * Polynomial poly2({1, 1}); // 1 + x
 * std::pair<Polynomial, Polynomial> result = poly1 / poly2;
 * Polynomial quotient = result.first;
 * Polynomial remainder = result.second;
 * // poly1 = poly2 * quotient + remainder
 * @endcode
 */
std::pair<Polynomial, Polynomial> Polynomial::operator/(const Polynomial &rhs) const
{
    return this->Divide(rhs);
}

/**
 * @brief Get all the roots of the Polynomial by Newton's method.
 * 
 * @return std::vector<double> 
 * 
 * @code {.cc}
 * Polynomial poly({1, 2, 1}); // 1 + 2x + x^2
 * std::vector<double> roots = poly.Get_all_roots(); // -1
 * @endcode
 */
std::vector<double> Polynomial::Get_all_roots() const
{
    std::vector<double> roots;
    Polynomial poly = *this;
    while (poly.get_degree() > 0)
    {
        NewtonSolver<Polynomial> newton(poly, 0);
        double root = newton.solve();
        if (std::abs(newton(root)) > 1e-3)
        {
            break;
        }
        else
        {
            roots.push_back(root);
            std::pair<Polynomial, Polynomial> result = poly.Divide(Polynomial({-root, 1}));
            poly = result.first;
        }
    }
    poly = *this;
    int n = roots.size();
    for (int i=0; i<n; i++)
    {
        NewtonSolver<Polynomial> newton(poly, roots[i], 1e-12, 1e-12, 1);
        roots[i] = newton.solve();
    }
    return roots;
}