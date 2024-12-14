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

/**
 * @brief This class is used to implement the curve interpolation, which can be
 * used to interpolate the curve in 2D, 3D, and the surface of a ball.
 * The curve is represented by a set of points, and the interpolation method is
 * determined by the method parameter.
 *
 * The curve can be interpolated by using the Bspline or PPoly method.
 *
 * @tparam N Order of the spline to be used.
 * @tparam Real The type of the data. Default is double.
 */
template <int N, typename Real = double>
class Curve
{
private:
    /**
     * @brief The knots of the curve
     */
    std::vector<Real> _t;
    /**
     * @brief The values of the curve
     */
    std::vector<std::vector<Real>> _xy;
    /**
     * @brief The method of the curve to be interpolated
     *
     * Related to the method of BInterpolate and PPInterpolate.
     */
    int _method;
    /**
     * @brief The boundary condition of the curve
     */
    std::vector<std::vector<Real>> _boundary_condition;

public:
    Curve() = default;
    /**
     * @brief Construct a new Curve object
     *
     * @param t The knots of the curve
     * @param xy The values of the curve
     * @param method The method of the curve to be interpolated
     * @param boundary_condition The boundary condition of the curve
     */
    Curve(const std::vector<Real>              &t,
          const std::vector<std::vector<Real>> &xy,
          const int                             method             = 0,
          const std::vector<std::vector<Real>> &boundary_condition = std::vector<std::vector<Real>>());
    /**
     * @brief Sqrt function
     * @param x The value to calculate the square root
     * @return Real The square root of x
     */
    Real
    SQRT(const Real &x) const; // the square root function
    /**
     * @brief Arcsin function
     * @param x The value to calculate the arcsin
     * @return Real The arcsin of x
     */
    Real
    ARCSIN(const Real &x) const; // the arcsin function
    /**
     * @brief Arctan function
     * @param x The value to calculate the arctan
     * @return Real The arctan of x
     */
    Real
    ARCCOS(const Real &x) const; // the arccos function
    /**
     * @brief Arctan function
     * @param x The value to calculate the arctan
     * @return Real The arctan of x
     */
    Real
    ARCTAN(const Real &x) const; // the arctan function
    /**
     * @brief Sin function
     * @param x The value to calculate the sin
     * @return Real The sin of x
     */
    Real
    SIN(const Real &x) const; // the sin function
    /**
     * @brief Cos function
     * @param x The value to calculate the cos
     * @return Real The cos of x
     */
    Real
    COS(const Real &x) const; // the cos function

    /**
     * @brief Use Bspline to interpolate the curve in 2D
     * @return std::vector<BSpline<Real>>
     */
    std::vector<BSpline<Real>>
    BSpline2d() const; // Use BSpline to interpolate
                       // the curve in 2D
    /**
     * @brief Use Bspline to interpolate the curve in 3D
     * @return std::vector<BSpline<Real>>
     */
    std::vector<BSpline<Real>>
    BSpline3d() const; // Use BSpline to interpolate
                       // the curve in 3D
    /**
     * @brief Use Bspline to interpolate the curve in the surface of a ball
     * @return BallFunction<Real>
     *
     * @details Write all points into
     * \f[
     * x = r \sin \theta \cos \phi, \quad
     * y = r \sin \theta \sin \phi, \quad
     * z = r \cos \theta
     * \f]
     *
     * Then use Bspline to interpolate \f$ \theta \f$ and \f$ \phi \f$.
     */
    BallFunction<Real>
    BSplineBall() const; // Use BSpline to interpolate
                         // the curve in the surface of a ball
    /**
     * @brief Use Bspline to interpolate the curve in the surface of a ball
     * @return BallFunction<Real>
     *
     * @details We can find a point \f$(x_0,y_0,z_0)\f$ not in the cruve, then the point can project the curve to a plane.
     * \f[
     * x_1 = \frac{z x_0 - z_0 x}{z - z_0}, \quad
     * y_1 = \frac{z y_0 - z_0 y}{z - z_0}
     * \f]
     *
     * Use Bspline to interpolate the projected curve.
     *
     * Finally, use the inverse formula to project the curve back to the surface of the ball.
     *
     * \f[
     * t = \frac{2 x_0 x + 2 y_0 y - 2 r^2}{r^2 - 2 x_0 x - 2 y_0 y + x^2 + y^2}, \quad
     * x_1 = x_0 + t (x_0 - x), \quad
     * y_1 = y_0 + t (y_0 - y), \quad
     * z_1 = z_0 + t (z_0 - z)
     * \f]
     *
     */
    BallFunction<Real>
    BSplineBallProj() const;

    /**
     * @brief Use PPoly to interpolate the curve in 2D
     * @return std::vector<PPoly<Real>>
     */
    std::vector<PPoly<Real>>
    PPoly2d() const; // Use PPoly to interpolate
                     // the curve in 2D

    /**
     * @brief Use PPoly to interpolate the curve in 3D
     * @return std::vector<PPoly<Real>>
     */
    std::vector<PPoly<Real>>
    PPoly3d() const; // Use PPoly to interpolate
                     // the curve in 3D
    /**
     * @brief Use PPoly to interpolate the curve in the surface of a ball
     * @return BallFunction<Real>
     *
     * @details Write all points into
     * \f[
     * x = r \sin \theta \cos \phi, \quad
     * y = r \sin \theta \sin \phi, \quad
     * z = r \cos \theta
     * \f]
     *
     * Then use PPoly to interpolate \f$ \theta \f$ and \f$ \phi \f$.
     */
    BallFunction<Real>
    PPolyBall() const; // Use PPoly to interpolate

    /**
     * @brief Use PPoly to interpolate the curve in the surface of a ball
     * @return BallFunction<Real>
     *
     * @details We can find a point \f$(x_0,y_0,z_0)\f$ not in the cruve, then the point can project the curve to a plane.
     * \f[
     * x_1 = \frac{z x_0 - z_0 x}{z - z_0}, \quad
     * y_1 = \frac{z y_0 - z_0 y}{z - z_0}
     * \f]
     *
     * Use PPoly to interpolate the projected curve.
     *
     * Finally, use the inverse formula to project the curve back to the surface of the ball.
     *
     * \f[
     * t = \frac{2 x_0 x + 2 y_0 y - 2 r^2}{r^2 - 2 x_0 x - 2 y_0 y + x^2 + y^2}, \quad
     * x_1 = x_0 + t (x_0 - x), \quad
     * y_1 = y_0 + t (y_0 - y), \quad
     * z_1 = z_0 + t (z_0 - z)
     * \f]
     *
     */
    BallFunction<Real>
    PPolyBallProj() const;

    /**
     * @brief Check if a point is in the curve
     * @param x x-coordinate of the point
     * @param y y-coordinate of the point
     * @param z z-coordinate of the point
     * @return true if the point is in the curve
     * @return false if the point is not in the curve
     */
    bool
    isPointInSet(Real x, Real y, Real z) const;

    /**
     * @brief Find a random point on the surface of the ball, in order to find a point not in the curve
     * @param x_0 x-coordinate of the point
     * @param y_0 y-coordinate of the point
     * @param z_0 z-coordinate of the point
     */
    void
    FindRandomPoint(Real &x_0, Real &y_0, Real &z_0) const;

    /**
     * @brief Calculate the theta and phi of the curve
     *
     * @param theta theta of the curve
     * @param phi phi of the curve
     */
    void
    CalculateThetaPhi(std::vector<Real> &theta, std::vector<Real> &phi) const;
};

#pragma once
#include "Curve.tpp"

#endif // CURVE_HPP