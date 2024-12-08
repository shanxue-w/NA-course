/**
 * @file Curve.tpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief This is the implementation of Curve
 * @version 0.1
 * @date 2024-11-24
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "BInterpolate.hpp"
#include "BSpline.hpp"
#include "Curve.hpp"

template <int N, typename Real>
Real
Curve<N, Real>::SQRT(const Real &x) const
{
    // if Real is double, float or long double, use the standard sqrt function
    // if Real is mpf_class, use the sqrt function in the gmp library

    if constexpr (
        std::is_same_v<
            Real,
            double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::sqrt(x);
    }
    else if constexpr (std::is_same_v<Real,
                                      mpf_class> || std::is_same_v<Real, mpf_t>)
    {
        Real result;
        mpf_sqrt(result.get_mpf_t(), x.get_mpf_t());
        return result;
    }
    else
    {
        return Real(0.0);
    }
}

template <int N, typename Real>
Real
Curve<N, Real>::ARCSIN(const Real &x) const
{
    if constexpr (
        std::is_same_v<
            Real,
            double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::asin(x);
    }
    else if constexpr (std::is_same_v<Real,
                                      mpf_class> || std::is_same_v<Real, mpf_t>)
    {
        return Real(std::asin(x.get_d()));
    }
    else
    {
        return Real(0.0);
    }
}

template <int N, typename Real>
Real
Curve<N, Real>::ARCCOS(const Real &x) const
{
    if constexpr (
        std::is_same_v<
            Real,
            double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::acos(x);
    }
    else if constexpr (std::is_same_v<Real,
                                      mpf_class> || std::is_same_v<Real, mpf_t>)
    {
        return Real(std::acos(x.get_d()));
    }
    else
    {
        return Real(0.0);
    }
}

template <int N, typename Real>
Real
Curve<N, Real>::ARCTAN(const Real &x) const
{
    if constexpr (
        std::is_same_v<
            Real,
            double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::atan(x);
    }
    else if constexpr (std::is_same_v<Real,
                                      mpf_class> || std::is_same_v<Real, mpf_t>)
    {
        return Real(std::atan(x.get_d()));
    }
    else
    {
        return Real(0.0);
    }
}

template <int N, typename Real>
Real
Curve<N, Real>::SIN(const Real &x) const
{
    if constexpr (
        std::is_same_v<
            Real,
            double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::sin(x);
    }
    else if constexpr (std::is_same_v<Real,
                                      mpf_class> || std::is_same_v<Real, mpf_t>)
    {
        return Real(std::sin(x.get_d()));
    }
    else
    {
        return Real(0.0);
    }
}

template <int N, typename Real>
Real
Curve<N, Real>::COS(const Real &x) const
{
    if constexpr (
        std::is_same_v<
            Real,
            double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::cos(x);
    }
    else if constexpr (std::is_same_v<Real,
                                      mpf_class> || std::is_same_v<Real, mpf_t>)
    {
        return Real(std::cos(x.get_d()));
    }
    else
    {
        return Real(0.0);
    }
}

template <int N, typename Real>
Curve<N, Real>::Curve(const std::vector<Real>              &t,
                      const std::vector<std::vector<Real>> &xy,
                      const int                             method,
                      const std::vector<std::vector<Real>> &boundary_condition)
    : _t(t), _xy(xy), _method(method)
{
    // if boundary_condition is not given, three dimensions are all set to 0
    if (boundary_condition.empty())
    {
        _boundary_condition = {{0}, {0}, {0}};
    }
}

template <int N, typename Real>
std::vector<BSpline<Real>>
Curve<N, Real>::BSpline2d() const
{
    BInterpolate<N, Real> b_x(_t, _xy[0], _method, _boundary_condition[0]); // x
    BInterpolate<N, Real> b_y(_t, _xy[1], _method, _boundary_condition[1]); // y
    BSpline<Real>         bs_x = b_x.getBSpline();
    BSpline<Real>         bs_y = b_y.getBSpline();
    return {bs_x, bs_y};
}

template <int N, typename Real>
std::vector<PPoly<Real>>
Curve<N, Real>::PPoly2d() const
{
    PPInterpolate<N, Real> p_x(
        _t, _xy[0], _method, _boundary_condition[0]); // x
    PPInterpolate<N, Real> p_y(
        _t, _xy[1], _method, _boundary_condition[1]); // y
    PPoly<Real> pp_x = p_x.getPoly();
    PPoly<Real> pp_y = p_y.getPoly();
    return {pp_x, pp_y};
}

template <int N, typename Real>
std::vector<BSpline<Real>>
Curve<N, Real>::BSpline3d() const
{
    BInterpolate<N, Real> b_x(_t, _xy[0], _method, _boundary_condition[0]); // x
    BInterpolate<N, Real> b_y(_t, _xy[1], _method, _boundary_condition[1]); // y
    BInterpolate<N, Real> b_z(_t, _xy[2], _method, _boundary_condition[2]); // z
    BSpline<Real>         bs_x = b_x.getBSpline();
    BSpline<Real>         bs_y = b_y.getBSpline();
    BSpline<Real>         bs_z = b_z.getBSpline();
    return {bs_x, bs_y, bs_z};
}

template <int N, typename Real>
std::vector<PPoly<Real>>
Curve<N, Real>::PPoly3d() const
{
    PPInterpolate<N, Real> p_x(
        _t, _xy[0], _method, _boundary_condition[0]); // x
    PPInterpolate<N, Real> p_y(
        _t, _xy[1], _method, _boundary_condition[1]); // y
    PPInterpolate<N, Real> p_z(
        _t, _xy[2], _method, _boundary_condition[2]); // z
    PPoly<Real> pp_x = p_x.getPoly();
    PPoly<Real> pp_y = p_y.getPoly();
    PPoly<Real> pp_z = p_z.getPoly();
    return {pp_x, pp_y, pp_z};
}

template <int N, typename Real>
BallFunction<Real>
Curve<N, Real>::BSplineBall() const
{
    /**
     * By polar coordinates, we can transform the 3D curve to a 2D curve
     *
     * \f$ x = r \sin \theta \cos \varphi \f$
     * \f$ y = r \sin \theta \sin \varphi \f$
     * \f$ z = r \cos \theta \f$
     */
    // initial radius
    Real r = SQRT(_xy[0][0] * _xy[0][0] + _xy[1][0] * _xy[1][0]
                  + _xy[2][0] * _xy[2][0]);
    // initial theta and phi
    int               t_size = _t.size();
    std::vector<Real> theta(t_size), phi(t_size);
    for (int i = 0; i < t_size; ++i)
    {
        theta[i]       = ARCCOS(_xy[2][i] / r);
        Real sin_theta = SIN(theta[i]);
        Real tmp;
        if (sin_theta > 1e-8)
        {
            // tmp = ARCCOS(_xy[0][i] / (r * sin_theta));
            // // if tmp is nan, it means that x = 0
            // if (std::isnan(tmp))
            // {
            //     tmp = 0;
            // }
            Real division = _xy[0][i] / (r * sin_theta);
            if (division > 1)
            {
                tmp = Real(0);
            }
            else if (division < -1)
            {
                tmp = Real(3.14159265358979323846) / Real(2);
            }
            else
            {
                tmp = ARCCOS(division);
            }
            if (_xy[1][i] < 0)
            {
                const Real PI = Real(3.14159265358979323846);
                tmp           = 2 * PI - tmp;
            }
            phi[i] = tmp;
        }
        else
        {
            phi[i] = 0;
        }
    }
    BSpline<Real> bs_theta =
        BInterpolate<N, Real>(_t, theta, _method, _boundary_condition[0])
            .getBSpline();
    BSpline<Real> bs_phi =
        BInterpolate<N, Real>(_t, phi, _method, _boundary_condition[0])
            .getBSpline();

    std::vector<BSpline<Real>> bs = {bs_theta, bs_phi};
    return BallFunction<Real>(bs, r);
}

template <int N, typename Real>
BallFunction<Real>
Curve<N, Real>::PPolyBall() const
{
    Real r = SQRT(_xy[0][0] * _xy[0][0] + _xy[1][0] * _xy[1][0]
                  + _xy[2][0] * _xy[2][0]);
    // initial theta and phi
    int               t_size = _t.size();
    std::vector<Real> theta(t_size), phi(t_size);
    for (int i = 0; i < t_size; ++i)
    {
        theta[i]       = ARCCOS(_xy[2][i] / r);
        Real sin_theta = SIN(theta[i]);
        Real tmp;
        if (sin_theta > 1e-8)
        {
            // tmp = ARCCOS(_xy[0][i] / (r * sin_theta));
            // // if tmp is nan, it means that x = 0
            // if (std::isnan(tmp))
            // {
            //     tmp = 0;
            // }
            Real division = _xy[0][i] / (r * sin_theta);
            if (division > 1)
            {
                tmp = Real(0);
            }
            else if (division < -1)
            {
                tmp = Real(3.14159265358979323846) / Real(2);
            }
            else
            {
                tmp = ARCCOS(division);
            }
            if (_xy[1][i] < 0)
            {
                const Real PI = Real(3.14159265358979323846);
                tmp           = 2 * PI - tmp;
            }
            phi[i] = tmp;
        }
        else
        {
            phi[i] = 0;
        }
    }
    PPoly<Real> pp_theta =
        PPInterpolate<N, Real>(_t, theta, _method, _boundary_condition[0])
            .getPoly();
    PPoly<Real> pp_phi =
        PPInterpolate<N, Real>(_t, phi, _method, _boundary_condition[0])
            .getPoly();

    std::vector<PPoly<Real>> pp = {pp_theta, pp_phi};
    return BallFunction<Real>(pp, r);
}