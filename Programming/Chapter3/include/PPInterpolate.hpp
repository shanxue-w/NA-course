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
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#define EIGEN_USE_THREADS

template <int N, typename Real = double> // N is the order of the polynomial and
                                         // Real is the type of the coefficients
class PPInterpolate {
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

  PPInterpolate(
      const std::vector<Real> &t, // nodes
      const std::vector<Real> &y, // values
      // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
      const int                method = 0,
      const std::vector<Real> &boundary_condition = std::vector<Real>(N, 0.0),
      const int                check = 0);

  void interpolate(
      const std::vector<Real> &t, // nodes
      const std::vector<Real> &y, // values
      // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
      const int                method = 0,
      const std::vector<Real> &boundary_condition = std::vector<Real>(N, 0.0));

  Real operator()(Real x) const;

  PPoly<Real> getPoly() const;
};

template <>
void PPInterpolate<1, double>::interpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot
    const std::vector<double> &boundary_condition); // boundary condition

template <>
void PPInterpolate<2, double>::interpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot
    const std::vector<double> &boundary_condition); // boundary condition

template <>
void PPInterpolate<3, double>::interpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot
    const std::vector<double> &boundary_condition); // boundary condition

#include "PPInterpolate.tpp"

#endif // PPINTERPOLATE_HPP