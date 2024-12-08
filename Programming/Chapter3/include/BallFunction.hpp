#ifndef __BALLFUNCTION_HPP__
#define __BALLFUNCTION_HPP__

#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"

template <typename Real>
class BallFunction
{
public:
    BallFunction(const std::vector<PPoly<Real>> &pp_ball_function,
                 const Real                     &r)
    {
        this->_pp_ball_function = pp_ball_function;
        this->r                 = r;
        this->method            = 0;
    }
    BallFunction(const std::vector<BSpline<Real>> &bs_ball_function,
                 const Real                       &r)
    {
        this->_bs_ball_function = bs_ball_function;
        this->r                 = r;
        this->method            = 1;
    }
    Real
    SIN(const Real &x) const
    {
        if constexpr (
            std::is_same_v<
                Real,
                double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
        {
            return std::sin(x);
        }
        else if constexpr (std::is_same_v<
                               Real,
                               mpf_class> || std::is_same_v<Real, mpf_t>)
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
        if constexpr (
            std::is_same_v<
                Real,
                double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
        {
            return std::cos(x);
        }
        else if constexpr (std::is_same_v<
                               Real,
                               mpf_class> || std::is_same_v<Real, mpf_t>)
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
            y[0] =
                r * SIN(_pp_ball_function[0](x)) * COS(_pp_ball_function[1](x));
            y[1] =
                r * SIN(_pp_ball_function[0](x)) * SIN(_pp_ball_function[1](x));
            y[2] = r * COS(_pp_ball_function[0](x));
        }
        else
        {
            y[0] =
                r * SIN(_bs_ball_function[0](x)) * COS(_bs_ball_function[1](x));
            y[1] =
                r * SIN(_bs_ball_function[0](x)) * SIN(_bs_ball_function[1](x));
            y[2] = r * COS(_bs_ball_function[0](x));
        }
        return y;
    }

private:
    std::vector<PPoly<Real>>   _pp_ball_function;
    std::vector<BSpline<Real>> _bs_ball_function;
    Real                       r;
    int                        method;
};

#endif