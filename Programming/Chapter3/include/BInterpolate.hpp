/**
 * @file BInterpolate.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Use B-spline to interpolate the data
 * @version 0.1
 * @date 2024-11-14
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef BINTERPOLATE_HPP
#define BINTERPOLATE_HPP

#include "BSpline.hpp"
#include <json/json.h>

#define EIGEN_USE_THREADS
#define EIGEN_USE_LAPACKE

template <int N, typename Real = double>
class BInterpolate
{
private:
    /**
     * nodes of the interpolation
     */
    std::vector<Real> _t;

    /**
     * interpolation data
     */
    std::vector<Real> _y;

    /**
     * The method of the interpolation
     */
    int _method;

    /**
     * Boundary condition for the interpolation, if needed.
     */
    std::vector<Real> _boundary_condition;

    /**
     * The result of the interpolation
     */
    BSpline<Real> _bspline;

public:
    BInterpolate() = default;

    BInterpolate(
        const std::vector<Real> &t,
        const std::vector<Real> &y,
        const int               &method             = 0,
        const std::vector<Real> &boundary_condition = std::vector<Real>(N, 0.0),
        const int                check              = 0);

    BInterpolate(const Json::Value &json); // load from json file

    void
    interpolate(const std::vector<Real> &t,
                const std::vector<Real> &y,
                const int               &method,
                const std::vector<Real> &boundary_condition =
                    std::vector<Real>(N, 0.0));

    BSpline<Real>
    getBSpline() const;

    Real
    operator()(const Real x);

    Real
    derivative(const Real x, const int n);
};

template <>
void
BInterpolate<1, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int                 &method,
    const std::vector<double> &boundary_condition);

template <>
void
BInterpolate<2, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int                 &method,
    const std::vector<double> &boundary_condition);

template <>
void
BInterpolate<3, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int &method, // 0 for periodic, 1 for complete, 2 for natural,
                       // 3 for not-a-knot.
    const std::vector<double> &boundary_condition);

template <>
void
BInterpolate<1, mpf_class>::interpolate(
    const std::vector<mpf_class> &t,
    const std::vector<mpf_class> &y,
    const int                    &method,
    const std::vector<mpf_class> &boundary_condition);

template <>
void
BInterpolate<2, mpf_class>::interpolate(
    const std::vector<mpf_class> &t,
    const std::vector<mpf_class> &y,
    const int                    &method,
    const std::vector<mpf_class> &boundary_condition);

template <>
void
BInterpolate<3, mpf_class>::interpolate(
    const std::vector<mpf_class> &t,
    const std::vector<mpf_class> &y,
    const int                    &method,
    const std::vector<mpf_class> &boundary_condition);

#include "BInterpolate.tpp"

#endif // BINTERPOLATE_HPP
