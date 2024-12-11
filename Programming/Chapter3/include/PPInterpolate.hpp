/**
 * @file PPInterpolate.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implmentation of PP-form interpolation
 * @version 0.1
 * @date 2024-11-09
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef PPINTERPOLATE_HPP
#define PPINTERPOLATE_HPP

#include "PPoly.hpp"
#include <json/json.h>
#include <type_traits>

#define EIGEN_USE_THREADS
#define EIGEN_USE_LAPACKE

/**
 * @brief Interpolation using `PPoly` class.
 *
 * @tparam N The order of the polynomial
 * @tparam Real The type of the coefficients
 *
 * @details The class PPInterpolate is used to interpolate the data using the
 * PP-form. The PP-form is a piecewise polynomial form, which is used to
 * represent a piecewise polynomial function. In each interval \f$ [t_{i},
 * t_{i+1}] \f$, the polynomial is represented as \f$ p_i(x) = \sum_{j=0}^{N-1}
 * a_{ij} (x-t_i)^j \f$.
 * Details about the ppform can be found in `PPoly.hpp`.
 *
 * The main thing this class does is to calculate the coefficients of the
 * polynomial in each interval.
 */
template <int N, typename Real = double> // N is the order of the polynomial and
                                         // Real is the type of the coefficients
class PPInterpolate
{
private:
    /**
     * nodes of the interpolation
     */
    std::vector<Real> _t;

    /**
     * \f$ y_i = f(t_i) \f$
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
    PPoly<Real> poly;

public:
    PPInterpolate() = default;

    PPInterpolate(const std::vector<Real> &t, // nodes
                  const std::vector<Real> &y, // values
                  // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
                  const int                method             = 0,
                  const std::vector<Real> &boundary_condition = std::vector<Real>(N, 0.0),
                  const int                check              = 0);

    PPInterpolate(const Json::Value &json); // load from json file

    void
    interpolate(const std::vector<Real> &t, // nodes
                const std::vector<Real> &y, // values
                // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
                const int                method             = 0,
                const std::vector<Real> &boundary_condition = std::vector<Real>(N, 0.0));

    Real
    operator()(Real x) const;

    PPoly<Real>
    getPoly() const;

protected:
    /**
     * @brief Generate the C matrix for the PP-form, which will be explained in
     * `interpolate` function.
     *
     * @return constexpr Eigen::Matrix<Real, N, N>
     */
    static constexpr Eigen::Matrix<Real, N, N>
    Generate_Cmatrix()
    {
        Eigen::Matrix<Real, N, N> C = Eigen::Matrix<Real, N, N>::Identity();
        for (int i = 1; i < N; i++)
        {
            C(0, i) = Real(i + 1);
        }
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                C(i, j) = C(i - 1, j - 1) + C(i, j - 1);
            }
        }
        return C;
    }
};

template <>
void
PPInterpolate<1, double>::interpolate(const std::vector<double> &t, // nodes
                                      const std::vector<double> &y, // values
                                      const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                        // not-a-knot
                                      const std::vector<double> &boundary_condition); // boundary condition

template <>
void
PPInterpolate<2, double>::interpolate(const std::vector<double> &t, // nodes
                                      const std::vector<double> &y, // values
                                      const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                        // not-a-knot
                                      const std::vector<double> &boundary_condition); // boundary condition

template <>
void
PPInterpolate<3, double>::interpolate(const std::vector<double> &t, // nodes
                                      const std::vector<double> &y, // values
                                      const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                        // not-a-knot
                                      const std::vector<double> &boundary_condition); // boundary condition

template <>
void
PPInterpolate<1, mpf_class>::interpolate(const std::vector<mpf_class> &t, // nodes
                                         const std::vector<mpf_class> &y, // values
                                         const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                           // not-a-knot
                                         const std::vector<mpf_class> &boundary_condition); // boundary condition

template <>
void
PPInterpolate<2, mpf_class>::interpolate(const std::vector<mpf_class> &t, // nodes
                                         const std::vector<mpf_class> &y, // values
                                         const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                           // not-a-knot
                                         const std::vector<mpf_class> &boundary_condition); // boundary condition

template <>
void
PPInterpolate<3, mpf_class>::interpolate(const std::vector<mpf_class> &t, // nodes
                                         const std::vector<mpf_class> &y, // values
                                         const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                                           // not-a-knot
                                         const std::vector<mpf_class> &boundary_condition); // boundary condition

#pragma once
#include "PPInterpolate.tpp"

#endif // PPINTERPOLATE_HPP