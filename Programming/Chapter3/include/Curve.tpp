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