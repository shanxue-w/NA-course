#include <iostream>
#include <vector>
#include "Hermite.hpp"
#include "NewtonPoly.hpp"
#include "Polynomial.hpp"


int main(void)
{
    std::vector<double> xData = {0, 3, 5, 8, 13, 0, 3, 5, 8, 13};
    std::vector<double> yData = {0, 225, 383, 623, 993, 75, 77, 80, 74, 72};
    std::vector<int> nData = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    Hermite hermite(xData, yData, nData);
    

    std::cout << "\n==================== D ====================" << std::endl;
    for (int i = 0; i < 10; i++)
    {
        if (i < 5)
        {
            std::cout << "Hermite(" << xData[i] << ") = " << hermite(xData[i]) << std::endl;
        }
        else
        {
            std::cout << "Hermite(" << xData[i] << ") = " << hermite.derivative(xData[i]) << std::endl;
        }
    }


    for (auto i=0.0; i<13.0; i+=0.01)
    {
        if (hermite.derivative(i) > 81.0)
        {
            std::cout << "Hermite'(" << i << ") = " << hermite.derivative(i) << std::endl;
            break;
        }
    }

    std::cout << "==================== D ====================\n" << std::endl;
    return 0;
}