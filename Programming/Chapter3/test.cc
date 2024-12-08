#include "PPInterpolate.hpp"
#include "PPoly.hpp"
#include <cmath>
#include <gmpxx.h>

int
main(void)
{
    int Max = 50;
    mpf_set_default_prec(3000);
    std::vector<mpf_class> x(Max);
    std::vector<mpf_class> y(Max);
    std::vector<mpf_class> boundary_conditions = {
        -1.0, 2.0, -6.0, 24.0, -120.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < Max; i++)
    {
        x[i] = i * i;
        y[i] = i;
    }

    PPInterpolate<4, mpf_class> inter(x, y, 0, boundary_conditions);
    PPoly<mpf_class>            poly = inter.getPoly();
    for (double i = 0; i < Max; i += 0.01)
    {
        std::cout << i << "," << inter(i) << std::endl;
    }
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

// #include "Curve.hpp"

// int
// main(void)
// {
// }