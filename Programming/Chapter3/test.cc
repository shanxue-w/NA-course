#include "PPInterpolate.hpp"
#include "PPoly.hpp"
#include <cmath>
#include <gmpxx.h>

int
main(void)
{
    double MAX = 100.0;
    // 设置精度为1000
    mpf_set_default_prec(10000);
    std::vector<mpf_class> x((int)MAX);
    std::vector<mpf_class> y((int)MAX);
    for (double i = 1.0; i <= MAX; i += 1.0)
    {
        x[i - 1] = i;
        y[i - 1] = i * i * i * i;
        // gmp log i
        // y[i-1] =
    }
    std::vector<mpf_class>      boundary = {4.0, 12.0, 24.0, 24.0};
    PPInterpolate<4, mpf_class> inter(x, y, 0, boundary, 0);
    PPoly<mpf_class>            poly = inter.getPoly();
    for (double i = 1.0; i <= MAX; i += 1)
    {
        std::cout << i << ", " << inter(i) << ", " << poly.derivative(i, 1)
                  << ", " << poly.derivative(i, 2) << ", "
                  << poly.derivative(i, 3) << std::endl;
    }
    return 0;
}

// #include "BSpline.hpp"

// int main(void) {
//   int                 MAX = 10;
//   std::vector<double> x(MAX);
//   std::vector<double> coeffs(MAX + 2, 0.0);
//   for (int i = 0; i < MAX; i++) {
//     x[i] = i;
//   }
//   for (int i = 0; i < MAX + 2; i += 1) {
//     coeffs[i] = 0.;
//   }
//   coeffs[5] = 1.0;
//   BSpline<double> bspline(coeffs, x, 3, 0);
//   for (double i = 0; i < 7.0; i += 1) {
//     // std::vector<double> basis = bspline.get_basis(i);
//     // for (auto j : basis) {
//     //   std::cout << j << ",";
//     // }
//     // std::cout << std::endl;
//     std::cout << i << "," << bspline(i) << std::endl;
//     // bspline(i);
//   }
// }

// #include "BInterpolate.hpp"
// #include "BSpline.hpp"
// #include <cmath>
// #include <gmpxx.h>

// int main(void) {
//   // mpf_set_default_prec(50);
//   // Eigen::initParallel();
//   // Eigen::setNbThreads(8);
//   int                 MAX = 300;
//   std::vector<double> x(MAX);
//   std::vector<double> y(MAX);
//   for (double i = 0; i < MAX; i++) {
//     x[i] = i;
//     y[i] = i + 1;
//   }
//   std::vector<double>      boundary = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//   0.0, 0.0,
//                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//                                        0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0,
//                                        0.0, 0.0, 0.0};
//   BInterpolate<10, double> inter(x, y, 0, boundary);
//   for (int i = 0; i < MAX; i++) {
//     std::cout << i << "," << inter(i) << "," << y[i] << std::endl;
//   }
//   for (int i = 1; i < 10; i++) {
//     std::cout << inter.derivative(0, i) << "," << inter.derivative(MAX - 1,
//     i)
//               << std::endl;
//   }
//   return 0;
// }