/**
 * @file Curve.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Use Bspline and PPoly to implement the curve interpolation
 * @version 0.1
 * @date 2024-11-22
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef CURVE_HPP
#define CURVE_HPP

#include "BInterpolate.hpp"
#include "BSpline.hpp"
#include "PPInterpolate.hpp"
#include <iostream>
#include <vector>

template <int N, typename Real = double> // N is the order of the curve and Real
                                         // is the type of the coefficients
class Curve {
private:
  std::vector<Real> _t; // the knots of the curve
  std::vector<Real> _x; // the x values of the curve
  std::vector<Real> _y; // the y values of the curve
  std::vector<Real> _z; // the z values of the curve, for N=3

public:
  Curve() = default;
  Curve(
      const std::vector<Real> &t,
      const std::vector<Real> &x,
      const std::vector<Real> &y,
      const std::vector<Real> &z = {}) :
      _t(t),
      _x(x), _y(y), _z(z) {}

  std::vector<BSpline<Real>> BSpline2d() const; // Use BSpline to interpolate
                                                // the curve in 2D
  std::vector<BSpline<Real>> BSpline3d() const; // Use BSpline to interpolate
                                                // the curve in 3D

  std::vector<PPoly<Real>> PPoly2d() const; // Use PPoly to interpolate
                                            // the curve in 2D
  std::vector<PPoly<Real>> PPoly3d() const; // Use PPoly to interpolate
                                            // the curve in 3D
};

#endif // CURVE_HPP