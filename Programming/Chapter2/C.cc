#include <iostream>
#include <vector>
#include "Newton.hpp"
#include "Lagrange.hpp"
#include <cmath>
#include "Polynomial.hpp"

double PI = M_PI;

std::vector<double> ChebyshevNodes(int n)
{
    std::vector<double> nodes = {};
    for (int i = 0; i < n; i++)
    {
        nodes.push_back(cos((2 * i + 1) * PI / (2 * n)));
    }
    return nodes;
}

double f(double x)
{
    return 1/(1+25.0*x*x);
}

int main(void)
{
    std::cout << "\n==================== C ====================" << std::endl;
    for (int n=5; n<=20; n+=5)
    {
        std::vector<double> xData = ChebyshevNodes(n);
        std::vector<double> yData = {};
        for (int i = 0; i < n; i++)
        {
            yData.push_back(f(xData[i]));
        }
        std::cout << "====================" << n << "====================" << std::endl;
        Newton newton(xData, yData);
        // std::cout << newton << std::endl;
        for (int i = 0; i < n; i++)
        {
            std::cout << "Newton(" << xData[i] << ") = " << newton(xData[i])
                      << " f(" << xData[i] << ") = " << f(xData[i]) << std::endl;
        }
        std::cout << "====================" << n << "====================" << std::endl;
    }
    std::cout << "==================== C ====================\ns" << std::endl;
}