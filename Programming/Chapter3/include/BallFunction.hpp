#ifndef __BALLFUNCTION_HPP__
#define __BALLFUNCTION_HPP__

#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"

template <typename Real>
class BallFunction
{
public:
    BallFunction(const std::vector<PPoly<Real>> &pp_ball_function, const Real &r)
    {
        this->_pp_ball_function = pp_ball_function;
        this->r                 = r;
        this->method            = 0;
    }
    BallFunction(const std::vector<BSpline<Real>> &bs_ball_function, const Real &r)
    {
        this->_bs_ball_function = bs_ball_function;
        this->r                 = r;
        this->method            = 1;
    }
    BallFunction(const std::vector<PPoly<Real>> &pp_ball_function, const std::vector<Real> &xyz, const Real &r)
    {
        this->_pp_ball_function = pp_ball_function;
        this->xyz               = xyz;
        this->r                 = r;
        this->method            = 2;
    }
    BallFunction(const std::vector<BSpline<Real>> &bs_ball_function, const std::vector<Real> &xyz, const Real &r)
    {
        this->_bs_ball_function = bs_ball_function;
        this->xyz               = xyz;
        this->r                 = r;
        this->method            = 3;
    }

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