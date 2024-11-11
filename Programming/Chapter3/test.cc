// #include "PPInterpolate.hpp"
// #include "PPoly.hpp"
// #include <cmath>


// int main(void)
// {
//     std::vector<double> x;
//     std::vector<double> y;
//     double MAX = 10000.0;
//     for (double i=1.0; i<=MAX; i+=1.0)
//     {
//         x.push_back(i);
//         y.push_back(std::log(i));
//     }
//     std::vector<double> boundary = {1.0, 1.0/MAX};
//     PPInterpolate<3> inter(x, y, 1, boundary, 0);
//     PPoly poly = inter.getPoly();
//     // for (double i=1.0; i<=MAX; i+=1.0)
//     // {
//     //     std::cout << "log(" << i << ") = " << std::log(i) << ", " << poly(i) << std::endl;
//     //     std::cout << "log'(" << i << ") = " << 1.0/i << ", " << poly.derivative(i, 1) << std::endl;
//     // }
//     return 0;
// }

#include "BSpline.hpp"


int main(void)
{
    int MAX = 1000;
    std::vector<double> x(MAX);
    std::vector<double> coeffs(MAX+9, 0.0);
    for (int i=0; i<MAX; i++)
    {
        x[i] = i;
    }
    coeffs[20] = 1.0;
    BSpline bspline(coeffs, x, 10);
    for (double i=9; i<21.0; i+=0.01)
    {
        std::cout << i << "," << bspline(i) << std::endl;
    }
}