/**
 * @file BallFunction.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implement the Ball Function Class
 * @version 0.1
 * @date 2024-12-15
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef __BALLFUNCTION_HPP__
#define __BALLFUNCTION_HPP__

#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"

/**
 * @brief Ball Function Class, which is used to calculate the ball function value.
 *
 * @tparam Real Type of the real number.
 */
template <typename Real>
class BallFunction
{
public:
    /**
     * @brief Construct a new Ball Function object.
     *
     * @param pp_ball_function The piecewise polynomial ball function.
     * @param r The radius of the ball.
     */
    BallFunction(const std::vector<PPoly<Real>> &pp_ball_function, const Real &r)
    {
        this->_pp_ball_function = pp_ball_function;
        this->r                 = r;
        this->method            = 0;
    }
    /**
     * @brief Construct a new Ball Function object
     *
     * @param bs_ball_function The B-spline ball function.
     * @param r The radius of the ball.
     */
    BallFunction(const std::vector<BSpline<Real>> &bs_ball_function, const Real &r)
    {
        this->_bs_ball_function = bs_ball_function;
        this->r                 = r;
        this->method            = 1;
    }
    /**
     * @brief Construct a new Ball Function object.
     *
     * @param pp_ball_function The piecewise polynomial ball function.
     * @param xyz Stereographic projection point.
     * @param r The radius of the ball.
     */
    BallFunction(const std::vector<PPoly<Real>> &pp_ball_function, const std::vector<Real> &xyz, const Real &r)
    {
        this->_pp_ball_function = pp_ball_function;
        this->xyz               = xyz;
        this->r                 = r;
        this->method            = 2;
    }
    /**
     * @brief Construct a new Ball Function object.
     *
     * @param bs_ball_function The B-spline ball function.
     * @param xyz Stereographic projection point.
     * @param r The radius of the ball.
     */
    BallFunction(const std::vector<BSpline<Real>> &bs_ball_function, const std::vector<Real> &xyz, const Real &r)
    {
        this->_bs_ball_function = bs_ball_function;
        this->xyz               = xyz;
        this->r                 = r;
        this->method            = 3;
    }

    /**
     * @brief Overload the sin function.
     *
     * @param x The input value.
     * @return Real The sin value.
     */
    Real
    SIN(const Real &x) const
    {
        if constexpr (std::is_same_v<Real, double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
        {
            return std::sin(x);
        }
        else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpf_t>)
        {
            return Real(std::sin(x.get_d()));
        }
        else
        {
            return Real(0.0);
        }
    }

    /**
     * @brief Overload the cos function.
     *
     * @param x The input value.
     * @return Real The cos value.
     */
    Real
    COS(const Real &x) const
    {
        if constexpr (std::is_same_v<Real, double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
        {
            return std::cos(x);
        }
        else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpf_t>)
        {
            return Real(std::cos(x.get_d()));
        }
        else
        {
            return Real(0.0);
        }
    }

    /**
     * @brief Calculate the ball function value.
     *
     * @param x The input value.
     * @return std::vector<Real> The ball function value.
     *
     * @details For Stereographic projection, we use the following formula:
     * \f[
     * t = \frac{2 x_0 x + 2 y_0 y - 2 r^2}{r^2 - 2 x_0 x - 2 y_0 y + x^2 + y^2}, \quad
     * x_1 = x_0 + t (x_0 - x), \quad
     * y_1 = y_0 + t (y_0 - y), \quad
     * z_1 = z_0 + t z_0.
     * \f]
     *
     * For the Spherical coordinates, we just write all points on the ball as
     * \f[
     * x = r \sin \theta \cos \phi, \quad
     * y = r \sin \theta \sin \phi, \quad
     * z = r \cos \theta.
     * \f]
     */
    std::vector<Real>
    operator()(const Real &x)
    {
        std::vector<Real> y(3);
        if (this->method == 0)
        {
            Real x_new = _pp_ball_function[0](x);
            Real y_new = _pp_ball_function[1](x);
            y[0]       = r * SIN(x_new) * COS(y_new);
            y[1]       = r * SIN(x_new) * SIN(y_new);
            y[2]       = r * COS(x_new);
        }
        else if (this->method == 1)
        {
            Real x_new = _bs_ball_function[0](x);
            Real y_new = _bs_ball_function[1](x);
            y[0]       = r * SIN(x_new) * COS(y_new);
            y[1]       = r * SIN(x_new) * SIN(y_new);
            y[2]       = r * COS(x_new);
        }
        else if (this->method == 2)
        {
            Real x_new = _pp_ball_function[0](x);
            Real y_new = _pp_ball_function[1](x);
            Real x0 = this->xyz[0], y0 = this->xyz[1], z0 = this->xyz[2];
            Real t = (2 * x0 * x_new + 2 * y0 * y_new - 2 * r * r)
                     / (r * r - 2 * x0 * x_new - 2 * y0 * y_new + x_new * x_new + y_new * y_new);
            y[0] = x0 + t * (x0 - x_new);
            y[1] = y0 + t * (y0 - y_new);
            y[2] = z0 + t * z0;
        }
        else if (this->method == 3)
        {
            Real x_new = _bs_ball_function[0](x);
            Real y_new = _bs_ball_function[1](x);
            Real x0 = this->xyz[0], y0 = this->xyz[1], z0 = this->xyz[2];
            Real t = (2 * x0 * x_new + 2 * y0 * y_new - 2 * r * r)
                     / (r * r - 2 * x0 * x_new - 2 * y0 * y_new + x_new * x_new + y_new * y_new);
            y[0] = x0 + t * (x0 - x_new);
            y[1] = y0 + t * (y0 - y_new);
            y[2] = z0 + t * z0;
        }
        return y;
    }

private:
    std::vector<PPoly<Real>>   _pp_ball_function;
    std::vector<BSpline<Real>> _bs_ball_function;
    Real                       r;
    int                        method;
    std::vector<Real>          xyz;
};

#endif