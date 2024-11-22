/**
 * @file BSpline.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Define the BSpline class.
 * @version 0.1
 * @date 2024-11-11
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef BSPLINE_HPP
#define BSPLINE_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <lapack.h>
#include <lapacke.h>
#include <omp.h>
#include <vector>

#define EIGEN_USE_LAPACKE
#define EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_THREADS

// using Real = double;

// template <typename T>
template <typename Real = double> // template <typename T>, default is double
class BSpline {
  /**
   * @brief The Bspline is in the form of
   * \f[ B(x) = \Sum_{j=i-n+1}^{i+1} c_j B_j^n(x), ~ x\in[x_i, x_{i+1}] \f]
   * So the computation of B_j^n(x) is very important.
   */
private:
  std::vector<Real>              _coeffs;
  std::vector<Real>              _t;
  std::vector<Real>              _basis;
  std::vector<std::vector<Real>> _total_basis;
  std::vector<Real>              _derivative_basis;
  int                            _n;

public:
  BSpline() = default;

  BSpline(
      const std::vector<Real> coeffs,
      const std::vector<Real> t,
      const int               n,
      const int               check = 0);

  BSpline(const BSpline &other);

  int get_interval(const Real x) const;

  BSpline &operator=(const BSpline &other);

  Real operator()(const Real x);

  std::vector<Real> get_coeffs() const;

  std::vector<Real> get_t() const;

  int get_n() const;

  std::vector<Real> get_basis(const Real x);

  void get_total_basis(const Real x);

  Real
  cal_derivative_basis(const int interval, const int j, const int degree, const int n);

  std::vector<Real> basis_derivative(const Real x, const int n);

  Real derivative(const Real x, const int n);
};

#include "BSpline.tpp"

#endif // BSPLINE_HPP