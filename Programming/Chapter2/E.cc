#include <iostream>
#include <vector>
#include "Newton.hpp"
#include "NewtonPoly.hpp"

int main(void)
{
    std::cout << "\n==================== E ====================" << std::endl;
    std::vector<double> xData = {0,6,10,13,17,20,28};
    std::vector<double> yData1= {6.67,17.3,42.7,37.3,30.1,29.3,28.7};
    std::vector<double> yData2= {6.67,16.1,18.9,15.0,10.6,9.44,8.89};

    Newton newton1(xData, yData1);
    Newton newton2(xData, yData2);

    Polynomial p1 = newton1.Convert_to_Polynomial();
    Polynomial p2 = newton2.Convert_to_Polynomial();

    std::vector<double> roots1 = p1.Get_all_roots();
    std::vector<double> roots2 = p2.Get_all_roots();
    
    for (auto i : xData)
    {
        std::cout << "Newton1(" << i << ") = " << newton1(i) << std::endl;
        std::cout << "Newton2(" << i << ") = " << newton2(i) << std::endl;
    }
    
    std::cout << "Roots of p1: " << std::endl;
    for (auto root : roots1)
    {
        std::cout << root << " ";
        std::cout << p1(root) << std::endl;
    }

    std::cout << "Roots of p2: ";
    for (auto root : roots2)
    {
        std::cout << root << " ";
    }
    std::cout << std::endl;



    std::cout << "==================== E ====================\n" << std::endl;

    return 0;
}