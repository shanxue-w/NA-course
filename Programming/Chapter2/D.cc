#include <iostream>
#include <vector>
#include "Hermite.hpp"
#include "Polynomial.hpp"

double newton(Polynomial poly, double x)
{
    while (std::abs(poly(x)) > 1e-12)
    {
        x = x - poly(x) / poly.derivative()(x);
    }
    return x;
}

int main()
{
    std::vector<double> xData = {0, 3, 5, 8, 13, 0, 3, 5, 8, 13};
    std::vector<double> yData = {0, 225, 383, 623, 993, 75, 77, 80, 74, 72};
    std::vector<int> nData = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    Hermite<Polynomial> hermite(xData, yData, nData);
    std::cout << hermite(10.0) << std::endl;

    Polynomial poly = hermite.get_polynomial(1);
    Polynomial poly2 = hermite.get_polynomial(2);
    std::vector<double> roots = poly2.Get_all_roots();
    for(auto root : roots)
    {
        std::cout << "t = " << root << ", speed = " << poly(root) << std::endl;
    }
    return 0;
}