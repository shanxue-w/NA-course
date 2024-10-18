#include "Polynomial.hpp"
#include "NewtonSolver.hpp"

Polynomial::Polynomial() : coefficients({0}), degree(0) {}

Polynomial::Polynomial(const std::vector<double> &coeffs) : coefficients(coeffs), degree(coeffs.size() - 1) {}

Polynomial::Polynomial(const Polynomial &rhs) : coefficients(rhs.coefficients), degree(rhs.degree) {}

Polynomial &Polynomial::operator=(const Polynomial &rhs)
{
    if (this != &rhs)
    {
        coefficients = rhs.coefficients;
        degree = rhs.degree;
    }
    return *this;
}

double Polynomial::operator()(double x) const
{
    double result = 0;
    for (int i = degree; i >= 0; i--)
    {
        result = result * x + coefficients[i];
    }
    return result;
}

double Polynomial::derivative(double x) const
{
    double result = 0;
    for (int i = degree; i >= 1; i--)
    {
        result = result * x + coefficients[i] * i;
    }
    return result;
}

double Polynomial::integral(double a, double b) const
{
    Polynomial integral = this->integral();
    return integral(b) - integral(a);
}

int Polynomial::get_degree() const
{
    return degree;
}

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

Polynomial Polynomial::operator*(double scalar) const
{
    std::vector<double> result(coefficients);
    for (int i = 0; i <= degree; i++)
    {
        result[i] *= scalar;
    }
    return Polynomial(result);
}

Polynomial Polynomial::operator/(double scalar) const
{
    std::vector<double> result(coefficients);
    for (int i = 0; i <= degree; i++)
    {
        result[i] /= scalar;
    }
    return Polynomial(result);
}

Polynomial Polynomial::derivative() const
{
    std::vector<double> result(std::max(degree - 1, 0) + 1, 0);
    for (int i = 1; i <= degree; i++)
    {
        result[i - 1] = coefficients[i] * i;
    }
    return Polynomial(result);
}

Polynomial Polynomial::integral() const
{
    std::vector<double> result(degree + 2, 0);
    for (int i = 0; i <= degree; i++)
    {
        result[i + 1] = coefficients[i] / (i + 1);
    }
    return Polynomial(result);
}

double Polynomial::integrate(double a, double b) const
{
    Polynomial integral = this->integral();
    return integral(b) - integral(a);
}

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

std::pair<Polynomial, Polynomial> Polynomial::operator/(const Polynomial &rhs) const
{
    return this->Divide(rhs);
}

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