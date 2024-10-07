#include <iostream>
#include <vector>
#include "Lagrange.hpp"


int main(void)
{
    std::vector<double> xData = {1,2,3,4};
    std::vector<double> yData = {1,4,7,11};

    Lagrange lagrange(xData, yData);
    std::cout << lagrange << std::endl;
    return 0;
}