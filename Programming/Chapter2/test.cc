#include <iostream>
#include <vector>
#include "Lagrange.hpp"
#include "Newton.hpp"


int main(void)
{
    std::vector<double> xData = {};
    std::vector<double> yData = {};
    for (int i = 0; i < 15; i++)
    {
        xData.push_back(i);
        yData.push_back(i * i);
    }

    Lagrange<Polynomial> lagrange(xData, yData);
    std::cout << lagrange << std::endl;

    Newton<Polynomial> newton(xData, yData);
    std::cout << newton << std::endl;
    return 0;
}