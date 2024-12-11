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

#include "BallFunction.hpp"

template <int N, typename Real = double> // N is the order of the curve and Real
                                         // is the type of the coefficients
class Curve
{
private:
    std::vector<Real>              _t;                  // the knots of the curve
    std::vector<std::vector<Real>> _xy;                 // the values of the curve
    int                            _method;             // the method of the curve
    std::vector<std::vector<Real>> _boundary_condition; // the boundary condition
public:
    Curve() = default;
    Curve(const std::vector<Real>              &t,
          const std::vector<std::vector<Real>> &xy,
          const int                             method             = 0,
          const std::vector<std::vector<Real>> &boundary_condition = std::vector<std::vector<Real>>());
    Real
    SQRT(const Real &x) const; // the square root function
    Real
    ARCSIN(const Real &x) const; // the arcsin function
    Real
    ARCCOS(const Real &x) const; // the arccos function
    Real
    ARCTAN(const Real &x) const; // the arctan function
    Real
    SIN(const Real &x) const; // the sin function
    Real
    COS(const Real &x) const; // the cos function

    std::vector<BSpline<Real>>
    BSpline2d() const; // Use BSpline to interpolate
                       // the curve in 2D
    std::vector<BSpline<Real>>
    BSpline3d() const; // Use BSpline to interpolate
                       // the curve in 3D
    BallFunction<Real>
    BSplineBall() const; // Use BSpline to interpolate
                         // the curve in the surface of a ball
    BallFunction<Real>
    BSplineBallProj() const;

    std::vector<PPoly<Real>>
    PPoly2d() const; // Use PPoly to interpolate
                     // the curve in 2D
    std::vector<PPoly<Real>>
    PPoly3d() const; // Use PPoly to interpolate
                     // the curve in 3D
    BallFunction<Real>
    PPolyBall() const; // Use PPoly to interpolate

    BallFunction<Real>
    PPolyBallProj() const;
};

#include "Curve.tpp"

#endif // CURVE_HPP