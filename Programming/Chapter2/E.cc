#include <iostream>
#include <vector>
#include "Newton.hpp"
#include "NewtonPoly.hpp"
#include <string>
#include <fstream>

int main(void)
{
    std::cout << "\n==================== E ====================" << std::endl;
    std::vector<double> xData = {0,6,10,13,17,20,28};
    std::vector<double> yData1= {6.67,17.3,42.7,37.3,30.1,29.3,28.7};
    std::vector<double> yData2= {6.67,16.1,18.9,15.0,10.6,9.44,8.89};

    Newton newton1(xData, yData1);
    Newton newton2(xData, yData2);

    Polynomial p1 = newton1.Convert_to_Polynomial().derivative();
    Polynomial p2 = newton2.Convert_to_Polynomial().derivative();

    std::vector<double> roots1 = p1.Get_all_roots();
    std::vector<double> roots2 = p2.Get_all_roots();
    

    
    std::cout << "Roots of p'1: " << std::endl;
    for (auto root : roots1)
    {
        std::cout << root << " ";
        std::cout << newton1(root) << std::endl;
    }

    std::cout << "Roots of p'2: " << std::endl;
    for (auto root : roots2)
    {
        std::cout << root << " ";
        std::cout << newton2(root) << std::endl;
    }
    std::cout << std::endl;


    std::cout << "==================== E ====================\n" << std::endl;

    std::string filename = "./data/E_Newton1.txt";
    std::ofstream file(filename);
    for (double j=0; j<=28+15; j+=0.01)
    {
        file << j << "," << newton1(j) << std::endl;
    }
    file.close();

    filename = "./data/E_Newton2.txt";
    std::ofstream file2(filename);
    for (double j=0; j<=28+15; j+=0.01)
    {
        file2 << j << "," << newton2(j) << std::endl;
    }
    file2.close();

    return 0;
}