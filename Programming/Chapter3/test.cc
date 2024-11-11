#include "PPInterpolate.hpp"
#include "PPoly.hpp"
#include <cmath>


int main(void)
{
    std::vector<double> x;
    std::vector<double> y;
    double MAX = 2.0;
    for (double i=1.0; i<=MAX; i+=1.0)
    {
        x.push_back(i);
        y.push_back(std::log(i));
    }
    std::vector<double> boundary = {1.0, 1.0/MAX};
    PPInterpolate<3> inter(x, y, 0, boundary);
    PPoly poly = inter.getPoly();
    for (double i=1.0; i<=MAX; i+=1.0)
    {
        std::cout << "log(" << i << ") = " << std::log(i) << ", " << poly(i) << std::endl;
        std::cout << "log'(" << i << ") = " << 1.0/i << ", " << poly.derivative(i, 1) << std::endl;
    }
    return 0;
}