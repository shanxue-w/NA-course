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
#include <lapack.h>
#include <lapacke.h>
#include <omp.h>
#include <string>
#include <vector>

#define EIGEN_USE_LAPACKE
#define EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_THREADS

template <int N> class PPInterpolate {
private:
  /**
   * nodes of the interpolation
   */
  std::vector<double> _t;

  /**
   * \f$ y_i = f(t_i) \f$
   */
  std::vector<double> _y;

  /**
   * The method of the interpolation
   */
  int _method;

  /**
   * Boundary condition for the interpolation, if needed.
   */
  std::vector<double> _boundary_condition;

  /**
   * The result of the interpolation
   */
  PPoly poly;

public:
  PPInterpolate() = default;

  PPInterpolate(
      const std::vector<double> &t,          // nodes
      const std::vector<double> &y,          // values
      const int                  method = 0, // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
      const std::vector<double> &boundary_condition = std::vector<double>(N, 0.0),
      const int                  check = 0);

  void interpolate(
      const std::vector<double> &t,          // nodes
      const std::vector<double> &y,          // values
      const int                  method = 0, // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
      const std::vector<double> &boundary_condition = std::vector<double>(N, 0.0));

  double operator()(double x) const;

  PPoly getPoly() const;
};

#endif // PPINTERPOLATE_HPP