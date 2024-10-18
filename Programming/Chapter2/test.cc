#include <iostream>
#include <vector>
#include "Hermite.hpp"
#include "Polynomial.hpp"

double f(double x)
{
    return 1 + x - 2*x*(x-1) + 2.0/3.0 *x*(x-1)*(x-1) - 5.0/36.0 *x*(x-1)*(x-1)*(x-3);
}

double f1(double x)
{
    return x*(x-1)*(x-1)*(x-3)*(x-3);
}

int main(void)
{
    // std::vector<double> xData = {0,1,3,1,3};
    // std::vector<double> yData = {1,2,0,-1,0};
    // std::vector<int> nData = {0,0,0,1,1};

    // Hermite<Polynomial> hermite(xData, yData, nData);
    // std::cout << hermite << std::endl;
    // std::cout << hermite(2) << std::endl;
    // std::cout << f(2) << std::endl;
    for (double i=0; i<=3; i+=0.01)
    {
        std::cout << i << " " << f1(i) << std::endl;
    }

    std::cout << (6+std::sqrt(21)) / 5.0 << " " << f1((6+std::sqrt(21)) / 5.0)/(1*2*3*4*5) << std::endl;
    return 0;
}