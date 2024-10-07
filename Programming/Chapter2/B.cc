#include <iostream>
#include <vector>
#include "Lagrange.hpp"
#include "Newton.hpp"
#include "Polynomial.hpp"

double f(double x)
{
    return 1/(1+x*x);
}

int main() 
{
    for (int n=2; n<=8; n+=2)
    {
        std::vector<double> xData = {};
        std::vector<double> yData = {};
        for (int i = 0; i < n; i++)
        {
            double l = -5 + 10 * (double)(i) / (double)(n);
            xData.push_back(l);
            yData.push_back(f(l));
        }
        std::cout << "====================" << n << "====================" << std::endl;
        Lagrange<Polynomial> lagrange(xData, yData);
        std::cout << lagrange << std::endl;

        Newton<Polynomial> newton(xData, yData);
        std::cout << newton << std::endl;
        std::cout << "====================" << n << "====================" << std::endl;
    }

    return 0;
}