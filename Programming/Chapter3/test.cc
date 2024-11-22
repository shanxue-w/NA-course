#include "PPInterpolate.hpp"
#include "PPoly.hpp"
#include <cmath>

int main(void) {
  double              MAX = 10000.0;
  std::vector<double> x((int)MAX);
  std::vector<double> y((int)MAX);
  for (double i = 1.0; i <= MAX; i += 1.0) {
    x[i - 1] = i;
    y[i - 1] = std::log(i);
  }
  std::vector<double>      boundary = {1.0, 1.0 / MAX};
  PPInterpolate<3, double> inter(x, y, 0, boundary, 1);
  PPoly<double>            poly = inter.getPoly();
  for (double i = 1.0; i <= MAX; i += 1.0) {
    std::cout << i << "," << inter(i) << "," << y[i - 1] << std::endl;
    std::cout << i << "," << poly.derivative(i) << "," << 1.0 / i << std::endl;
  }
  return 0;
}

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
//     for (int i=0; i<MAX+2; i+=1)
//     {
//         coeffs[i] = i;
//     }
//     coeffs[5] = 1.0;
//     BSpline bspline(coeffs, x, 3, 0);
//     for (double i=0; i<7.0; i+=1)
//     {
//         // std::vector<double> basis = bspline.get_basis(i);
//         // for (auto j: basis)
//         // {
//         //     std::cout << j << ",";
//         // }
//         // std::cout << std::endl;
//         // std::cout << i << "," << bspline(i) << std::endl;
//         bspline(i);
//     }
// }

// #include "BInterpolate.hpp"
// #include "BSpline.hpp"
// #include <gmpxx.h>

// int main(void) {
//   // mpf_set_default_prec(50);
//   // Eigen::initParallel();
//   // Eigen::setNbThreads(8);
//   int                 MAX = 3000;
//   std::vector<double> x(MAX);
//   std::vector<double> y(MAX);
//   for (double i = 0; i < MAX; i++) {
//     x[i] = i;
//     y[i] = i;
//   }
//   std::vector<double>      boundary = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//   0.0, 0.0,
//                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//                                        0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0,
//                                        0.0, 0.0, 0.0};
//   BInterpolate<10, double> inter(x, y, 1, boundary);
//   for (int i = 0; i < MAX; i++) {
//     std::cout << i << "," << inter(i) << "," << y[i] << std::endl;
//   }
//   return 0;
// }