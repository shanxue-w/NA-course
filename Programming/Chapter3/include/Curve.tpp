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

#ifndef CURVE_TPP
#define CURVE_TPP

#include "BInterpolate.hpp"
#include "BSpline.hpp"
#include "Curve.hpp"

template <int N, typename Real>
Real
Curve<N, Real>::SQRT(const Real &x) const
{
    // if Real is double, float or long double, use the standard sqrt function
    // if Real is mpf_class, use the sqrt function in the gmp library

    if constexpr (std::is_same_v<Real, double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::sqrt(x);
    }
    else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpf_t>)
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
    if constexpr (std::is_same_v<Real, double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::asin(x);
    }
    else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpf_t>)
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
    if constexpr (std::is_same_v<Real, double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::acos(x);
    }
    else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpf_t>)
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
    if constexpr (std::is_same_v<Real, double> or std::is_same_v<Real, float> or std::is_same_v<Real, long double>)
    {
        return std::atan(x);
    }
    else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpf_t>)
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

template <int N, typename Real>
Real
Curve<N, Real>::COS(const Real &x) const
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
    PPInterpolate<N, Real> p_x(_t, _xy[0], _method, _boundary_condition[0]); // x
    PPInterpolate<N, Real> p_y(_t, _xy[1], _method, _boundary_condition[1]); // y
    PPoly<Real>            pp_x = p_x.getPoly();
    PPoly<Real>            pp_y = p_y.getPoly();
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
    PPInterpolate<N, Real> p_x(_t, _xy[0], _method, _boundary_condition[0]); // x
    PPInterpolate<N, Real> p_y(_t, _xy[1], _method, _boundary_condition[1]); // y
    PPInterpolate<N, Real> p_z(_t, _xy[2], _method, _boundary_condition[2]); // z
    PPoly<Real>            pp_x = p_x.getPoly();
    PPoly<Real>            pp_y = p_y.getPoly();
    PPoly<Real>            pp_z = p_z.getPoly();
    return {pp_x, pp_y, pp_z};
}

// Read one of `BSplineBall`, `PPolyBall` is ok, just the same code.
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
    Real r = SQRT(_xy[0][0] * _xy[0][0] + _xy[1][0] * _xy[1][0] + _xy[2][0] * _xy[2][0]);
    // initial theta and phi
    int               t_size = _t.size();
    std::vector<Real> theta(t_size, Real(0)), phi(t_size, Real(0));
    int               flag = 0;
    for (int i = 0; i < t_size; ++i)
    {
        theta[i]       = ARCCOS(_xy[2][i] / r);
        Real sin_theta = SIN(theta[i]);
        Real tmp;
        if (sin_theta > 1e-10)
        {
            const Real PI       = Real(M_PI);
            Real       division = _xy[1][i] / _xy[0][i];
            tmp                 = ARCTAN(division);
            if (_xy[0][i] < 0)
            {
                tmp += PI;
            }
            phi[i] = tmp;
            if (flag == 1)
            {
                phi[i - 1] = phi[i];
                flag       = 0;
            }
            if (i > 0)
            {
                if (phi[i] - phi[i - 1] > 1)
                {
                    theta[i] = -theta[i];
                    phi[i]   = phi[i] - PI;
                }
                else if (phi[i] - phi[i - 1] < -1)
                {
                    theta[i] = -theta[i];
                    phi[i]   = phi[i] + PI;
                }
            }
        }
        else
        {
            if (i != 0)
            {
                phi[i] = phi[i - 1];
            }
            else
            {
                flag = 1;
            }
        }
    }
    BSpline<Real> bs_theta = BInterpolate<N, Real>(_t, theta, _method, _boundary_condition[0]).getBSpline();
    BSpline<Real> bs_phi   = BInterpolate<N, Real>(_t, phi, _method, _boundary_condition[0]).getBSpline();

    std::vector<BSpline<Real>> bs = {bs_theta, bs_phi};
    return BallFunction<Real>(bs, r);
}

// Read one of `BSplineBall`, `PPolyBall` is ok, just the same code.
template <int N, typename Real>
BallFunction<Real>
Curve<N, Real>::PPolyBall() const
{
    Real r = SQRT(_xy[0][0] * _xy[0][0] + _xy[1][0] * _xy[1][0] + _xy[2][0] * _xy[2][0]);
    // initial theta and phi
    int               t_size = _t.size();
    std::vector<Real> theta(t_size, Real(0)), phi(t_size, Real(0));
    int               flag = 0;
    for (int i = 0; i < t_size; ++i)
    {
        theta[i]       = ARCCOS(_xy[2][i] / r);
        Real sin_theta = SIN(theta[i]);
        Real tmp;
        if (sin_theta > 1e-10)
        {
            const Real PI       = Real(M_PI);
            Real       division = _xy[1][i] / _xy[0][i];
            tmp                 = ARCTAN(division);
            if (_xy[0][i] < 0)
            {
                tmp += PI;
            }
            phi[i] = tmp;
            if (flag == 1)
            {
                phi[i - 1] = phi[i];
                flag       = 0;
            }
            if (i > 0)
            {
                if (phi[i] - phi[i - 1] > 1)
                {
                    theta[i] = -theta[i];
                    phi[i]   = phi[i] - PI;
                }
                else if (phi[i] - phi[i - 1] < -1)
                {
                    theta[i] = -theta[i];
                    phi[i]   = phi[i] + PI;
                }
            }
        }
        else
        {
            if (i != 0)
            {
                phi[i] = phi[i - 1];
            }
            else
            {
                flag = 1;
            }
        }
    }
    PPoly<Real> pp_theta = PPInterpolate<N, Real>(_t, theta, _method, _boundary_condition[0]).getPoly();
    PPoly<Real> pp_phi   = PPInterpolate<N, Real>(_t, phi, _method, _boundary_condition[0]).getPoly();

    std::vector<PPoly<Real>> pp = {pp_theta, pp_phi};
    return BallFunction<Real>(pp, r);
}

// Read one of `BSplineBallProj`, `PPolyBallProj` is ok, just the same code.
template <int N, typename Real>
BallFunction<Real>
Curve<N, Real>::BSplineBallProj() const
{
    Real r = SQRT(_xy[0][0] * _xy[0][0] + _xy[1][0] * _xy[1][0] + _xy[2][0] * _xy[2][0]);
    Real x_0, y_0, z_0;
    // define the coordinates of the south pole and the north pole
    Real southPole[3] = {0.0, 0.0, -1.0};
    Real northPole[3] = {0.0, 0.0, 1.0};

    // auxiliary function to check if a point is in the set
    auto isPointInSet = [&](Real x, Real y, Real z)
    {
        if constexpr (std::is_same_v<Real, double> || std::is_same_v<Real, float> || std::is_same_v<Real, long double>)
        {
            for (int i = 0; i < N; ++i)
            {
                if (std::abs(_xy[0][i] - x) < 1e-6 && std::abs(_xy[1][i] - y) < 1e-6 && std::abs(_xy[2][i] - z) < 1e-6)
                {
                    return true;
                }
            }
        }
        else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpq_class>)
        {
            for (int i = 0; i < N; ++i)
            {
                if (abs(_xy[0][i] - x) < 1e-6 && abs(_xy[1][i] - y) < 1e-6 && abs(_xy[2][i] - z) < 1e-6)
                {
                    return true;
                }
            }
        }
        return false;
    };

    // check if the south pole and the north pole are occupied
    bool southOccupied = isPointInSet(southPole[0], southPole[1], southPole[2]);
    bool northOccupied = isPointInSet(northPole[0], northPole[1], northPole[2]);

    // if the south pole is not occupied, return the south pole
    if (!southOccupied)
    {
        x_0 = southPole[0];
        y_0 = southPole[1];
        z_0 = southPole[2];
    }
    // if the north pole is not occupied, return the north pole
    else if (!northOccupied)
    {
        x_0 = northPole[0];
        y_0 = northPole[1];
        z_0 = northPole[2];
    }
    else
    {
        // if both the south pole and the north pole are occupied, randomly select a point not in the set
        // here we generate a random point on the sphere in spherical coordinates
        bool found = false;
        while (!found)
        {
            // randomly generate a point on the sphere
            Real theta = 2.0 * M_PI * (rand() / Real(RAND_MAX));        // 随机生成 [0, 2π] 之间的角度
            Real phi   = ARCCOS(2.0 * (rand() / Real(RAND_MAX)) - 1.0); // 随机生成 [0, π] 之间的角度

            // convert to Cartesian coordinates
            Real x = r * SIN(phi) * COS(theta);
            Real y = r * SIN(phi) * SIN(theta);
            Real z = r * COS(phi);

            // check if the generated point already exists
            if (!isPointInSet(x, y, z))
            {
                x_0   = x;
                y_0   = y;
                z_0   = z;
                found = true;
            }
        }
    }
    int t_size = _t.size();

    std::vector<Real> x_proj(t_size), y_proj(t_size);
    for (int i = 0; i < t_size; ++i)
    {
        x_proj[i] = (_xy[2][i] * x_0 - _xy[0][i] * z_0) / (_xy[2][i] - z_0);
        y_proj[i] = (_xy[2][i] * y_0 - _xy[1][i] * z_0) / (_xy[2][i] - z_0);
    }

    BSpline<Real> bs_x = BInterpolate<N, Real>(_t, x_proj, _method, _boundary_condition[0]).getBSpline();
    BSpline<Real> bs_y = BInterpolate<N, Real>(_t, y_proj, _method, _boundary_condition[0]).getBSpline();

    std::vector<BSpline<Real>> bs     = {bs_x, bs_y};
    std::vector<Real>          center = {x_0, y_0, z_0};
    return BallFunction<Real>(bs, center, r);
}

// Read one of `BSplineBallProj`, `PPolyBallProj` is ok, just the same code.
template <int N, typename Real>
BallFunction<Real>
Curve<N, Real>::PPolyBallProj() const
{
    Real r = SQRT(_xy[0][0] * _xy[0][0] + _xy[1][0] * _xy[1][0] + _xy[2][0] * _xy[2][0]);
    Real x_0, y_0, z_0;
    // define the coordinates of the south pole and the north pole
    Real southPole[3] = {0.0, 0.0, -1.0};
    Real northPole[3] = {0.0, 0.0, 1.0};

    // auxiliary function to check if a point is in the set
    auto isPointInSet = [&](Real x, Real y, Real z)
    {
        if constexpr (std::is_same_v<Real, double> || std::is_same_v<Real, float> || std::is_same_v<Real, long double>)
        {
            for (int i = 0; i < N; ++i)
            {
                if (std::abs(_xy[0][i] - x) < 1e-6 && std::abs(_xy[1][i] - y) < 1e-6 && std::abs(_xy[2][i] - z) < 1e-6)
                {
                    return true;
                }
            }
        }
        else if constexpr (std::is_same_v<Real, mpf_class> || std::is_same_v<Real, mpq_class>)
        {
            for (int i = 0; i < N; ++i)
            {
                if (abs(_xy[0][i] - x) < 1e-6 && abs(_xy[1][i] - y) < 1e-6 && abs(_xy[2][i] - z) < 1e-6)
                {
                    return true;
                }
            }
        }
        return false;
    };

    // check if the south pole and the north pole are occupied
    bool southOccupied = isPointInSet(southPole[0], southPole[1], southPole[2]);
    bool northOccupied = isPointInSet(northPole[0], northPole[1], northPole[2]);

    // if the south pole is not occupied, return the south pole
    if (!southOccupied)
    {
        x_0 = southPole[0];
        y_0 = southPole[1];
        z_0 = southPole[2];
    }
    // if the north pole is not occupied, return the north pole
    else if (!northOccupied)
    {
        x_0 = northPole[0];
        y_0 = northPole[1];
        z_0 = northPole[2];
    }
    else
    {
        // if both the south pole and the north pole are occupied, randomly select a point not in the set
        // here we generate a random point on the sphere in spherical coordinates
        bool found = false;
        while (!found)
        {
            // randomly generate a point on the sphere
            Real theta = 2.0 * M_PI * (rand() / Real(RAND_MAX));        // 随机生成 [0, 2π] 之间的角度
            Real phi   = ARCCOS(2.0 * (rand() / Real(RAND_MAX)) - 1.0); // 随机生成 [0, π] 之间的角度

            // convert to Cartesian coordinates
            Real x = r * SIN(phi) * COS(theta);
            Real y = r * SIN(phi) * SIN(theta);
            Real z = r * COS(phi);

            // check if the generated point already exists
            if (!isPointInSet(x, y, z))
            {
                x_0   = x;
                y_0   = y;
                z_0   = z;
                found = true;
            }
        }
    }
    int t_size = _t.size();

    std::vector<Real> x_proj(t_size), y_proj(t_size);
    for (int i = 0; i < t_size; ++i)
    {
        x_proj[i] = (_xy[2][i] * x_0 - _xy[0][i] * z_0) / (_xy[2][i] - z_0);
        y_proj[i] = (_xy[2][i] * y_0 - _xy[1][i] * z_0) / (_xy[2][i] - z_0);
    }

    PPoly<Real> pp_x = PPInterpolate<N, Real>(_t, x_proj, _method, _boundary_condition[0]).getPoly();
    PPoly<Real> pp_y = PPInterpolate<N, Real>(_t, y_proj, _method, _boundary_condition[0]).getPoly();

    std::vector<PPoly<Real>> pp     = {pp_x, pp_y};
    std::vector<Real>        center = {x_0, y_0, z_0};
    return BallFunction<Real>(pp, center, r);
}

#endif // CURVE_TPP