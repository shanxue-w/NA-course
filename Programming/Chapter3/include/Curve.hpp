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

template <int N, typename Real = double> // N is the order of the curve and Real
                                         // is the type of the coefficients
class Curve
{
private:
    std::vector<Real>              _t;      // the knots of the curve
    std::vector<std::vector<Real>> _xy;     // the values of the curve
    int                            _method; // the method of the curve
    std::vector<std::vector<Real>>
        _boundary_condition; // the boundary condition
public:
    Curve() = default;
    Curve(const std::vector<Real>              &t,
          const std::vector<std::vector<Real>> &xy,
          const int                             method = 0,
          const std::vector<std::vector<Real>> &boundary_condition =
              std::vector<std::vector<Real>>());

    std::vector<BSpline<Real>>
    BSpline2d() const; // Use BSpline to interpolate
                       // the curve in 2D
    std::vector<BSpline<Real>>
    BSpline3d() const; // Use BSpline to interpolate
                       // the curve in 3D

    std::vector<PPoly<Real>>
    PPoly2d() const; // Use PPoly to interpolate
                     // the curve in 2D
    std::vector<PPoly<Real>>
    PPoly3d() const; // Use PPoly to interpolate
                     // the curve in 3D
};

#include "Curve.tpp"

#endif // CURVE_HPP