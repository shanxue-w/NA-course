// #include "PPInterpolate.hpp"
// #include "PPoly.hpp"
// #include <cmath>


// int main(void)
// {
//     double MAX = 10000.0;
//     std::vector<double> x((int)MAX);
//     std::vector<double> y((int)MAX);
//     for (double i=1.0; i<=MAX; i+=1.0)
//     {
//         x[i-1] = i;
//         y[i-1] = std::log(i);
//     }
//     std::vector<double> boundary = {1.0, 1.0/MAX};
//     PPInterpolate<3> inter(x, y, 0, boundary, 1);
//     PPoly poly = inter.getPoly();
//     return 0;
// }

// #include "BSpline.hpp"


// int main(void)
// {
//     int MAX = 10;
//     std::vector<double> x(MAX);
//     std::vector<double> coeffs(MAX+2, 0.0);
//     for (int i=0; i<MAX; i++)
//     {
//         x[i] = i;
//     }
//     coeffs[5] = 1.0;
//     BSpline bspline(coeffs, x, 3, 0);
//     for (double i=2; i<7.0; i+=0.01)
//     {
//         std::cout << i << "," << bspline.derivative(i, 4) << std::endl;
//         // std::cout << i << "," << bspline(i) << std::endl;
//     }
// }

#include "BInterpolate.hpp"
#include "BSpline.hpp"


int main(void)
{
    int MAX = 10;
    std::vector<double> x(MAX);
    std::vector<double> y(MAX);
    for (int i=0; i<MAX; i++)
    {
        x[i] = i;
        y[i] = std::sin(i);
    }
    BInterpolate<3> inter(x, y, 2);
    std::cout << inter.derivative(x[0], 1) << ", " << inter.derivative(x[MAX-1], 1) << std::endl;
    std::cout << inter.derivative(x[0], 2) << ", " << inter.derivative(x[MAX-1], 2) << std::endl;
    return 0;
}