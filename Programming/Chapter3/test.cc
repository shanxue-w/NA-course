#include "PPInterpolate.hpp"
#include "PPoly.hpp"
#include <cmath>


int main(void)
{
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 6.0};
    std::vector<double> y = {std::log(1.0), std::log(2.0), std::log(3.0), std::log(4.0), std::log(6.0)};
    std::vector<double> boundary = {1.0, 1.0/6.0};
    PPInterpolate<3> inter(x, y, 1, boundary);
    for (double i=1.0; i<=6.0; i+=0.01)
    {
        std::cout << i << " " << inter(i) << std::endl;
    }
}