#include <iostream>
#include <vector>
#include "Hermite.hpp"
#include "NewtonPoly.hpp"
#include "Polynomial.hpp"
#include <fstream>
#include <string>


int main(void)
{
    std::vector<double> xData = {0, 3, 5, 8, 13, 0, 3, 5, 8, 13};
    std::vector<double> yData = {0, 225, 383, 623, 993, 75, 77, 80, 74, 72};
    std::vector<int> nData = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    Hermite hermite(xData, yData, nData);
    

    std::cout << "\n==================== D ====================" << std::endl;
    std::cout << "Hermite(" << 10.0 << ") = " << hermite(10.0) << std::endl;

    std::cout << "\n============= Search by roots =============\n";
    Polynomial poly = hermite.Convert_to_Polynomial(); // Convert to polynomial
    Polynomial poly_prime = poly.derivative(); // Get the derivative of the polynomial
    Polynomial poly_2_prime = poly_prime.derivative(); // Get the second derivative of the polynomial
    std::vector<double> roots = poly_2_prime.Get_all_roots(); // Get all the roots of the polynomial
    for (auto root : roots)
    {
        std::cout << "Hermite'(" << root << ") = " << hermite.derivative(root) << std::endl;
    }

    std::cout << "\n============= Search by iteration =============\n";
    for (auto i=0.0; i<13.0; i+=0.01)
    {
        if (hermite.derivative(i) > 81.0)
        {
            std::cout << "Hermite'(" << i << ") = " << hermite.derivative(i) << std::endl;
            break;
        }
    }

    std::string filename = "./data/D_Hermite.txt";
    std::ofstream file(filename);
    for (double j=0; j<=13; j+=0.01)
    {
        file << j << "," << hermite(j) << std::endl;
    }
    file.close();

    std::string filename1 = "./data/D_Hermite_prime.txt";
    std::ofstream file1(filename1);
    for (double j=0; j<=13; j+=0.01)
    {
        file1 << j << "," << hermite.derivative(j) << std::endl;
    }
    file1.close();    

    std::cout << "==================== D ====================\n" << std::endl;
    return 0;
}