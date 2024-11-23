/**
 * @file PPInterpolate.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implemetnation of PPInterpolate class
 * @version 0.1
 * @date 2024-11-09
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "PPInterpolate.hpp"
#include <Eigen/src/Core/Matrix.h>

template <>
PPInterpolate<1, double>::PPInterpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot.
    const std::vector<double> &boundary_condition,
    const int                  check) // whether to check the order of t
    :
    _t(t),
    _y(y), _method(method), _boundary_condition(boundary_condition) {
  poly = PPoly<double>();
  if (check) {
    if (!std::is_sorted(t.begin(), t.end())) {
      std::vector<int> idx(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        idx[i] = i;
      }
      std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });
      // std::vector<Real> t_sorted(t.size()), y_sorted(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        _t[i] = t[idx[i]];
        _y[i] = y[idx[i]];
      }
    }
  }
  interpolate(t, y, method, boundary_condition); // interpolate
}

template <>
void PPInterpolate<1, double>::interpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot.
    const std::vector<double> &boundary_condition) // boundary condition
{
  // unused variable, how to solve it.
  (void)method;
  (void)boundary_condition;
  std::vector<std::vector<double>> A(t.size() - 1, std::vector<double>(2, 0));

  for (size_t i = 0; i < t.size() - 1; ++i) {
    A[i][0] = y[i];
    A[i][1] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
  }
  poly = PPoly<double>(A, t, 0); // interpolate
}

template <>
PPInterpolate<2, double>::PPInterpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot.
    const std::vector<double> &boundary_condition,
    const int                  check) // whether to check the order of t
    :
    _t(t),
    _y(y), _method(method), _boundary_condition(boundary_condition) {
  poly = PPoly<double>();
  if (check) {
    if (!std::is_sorted(t.begin(), t.end())) {
      std::vector<int> idx(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        idx[i] = i;
      }
      std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });
      // std::vector<Real> t_sorted(t.size()), y_sorted(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        _t[i] = t[idx[i]];
        _y[i] = y[idx[i]];
      }
    }
  }
  interpolate(t, y, method, boundary_condition); // interpolate
}

template <>
void PPInterpolate<2, double>::interpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot.
    const std::vector<double> &boundary_condition) // boundary condition
{
  switch (method) {
  case 0: {
    int                 t_size = t.size();
    std::vector<double> K_i(t_size - 1, 0);

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    Eigen::SparseMatrix<double> A(t_size, t_size);
    Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);

    std::vector<Eigen::Triplet<double>>
        triplets; // For batch insertion into the sparse matrix

    for (int i = 0; i < t_size - 1; i++) {
      b(i) = 2 * K_i[i];

      {
        triplets.push_back(Eigen::Triplet<double>(i, i, 1));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, 1));
      }
    }
    triplets.push_back(Eigen::Triplet<double>(t_size - 1, 0, 1));
    triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 1, -1));
    b(t_size - 1) = 0;
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXd c = solver.solve(b);

    std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(3, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      coeffs[i][0] = y[i];
      coeffs[i][1] = c(i);
      coeffs[i][2] = (c(i + 1) - c(i)) / (2 * (t[i + 1] - t[i]));
    }
    poly = PPoly<double>(coeffs, t, 0); // interpolate
    return;
  } break;

  case 1: {
    int                 t_size = t.size();
    std::vector<double> K_i(t_size - 1, 0);

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    std::vector<double> c(t_size, 0);
    c[0] = boundary_condition[0];

    for (int i = 1; i < t_size; ++i) {
      c[i] = 2 * K_i[i - 1] - c[i - 1];
    }

    std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(3, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      coeffs[i][0] = y[i];
      coeffs[i][1] = c[i];
      coeffs[i][2] = (c[i + 1] - c[i]) / (2 * (t[i + 1] - t[i]));
    }
    poly = PPoly<double>(coeffs, t, 0); // interpolate
    return;
  } break;

  default:
    throw std::invalid_argument("method must be 0,1");
    return;
    break;
  }
  return;
}

template <>
PPInterpolate<3, double>::PPInterpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot.
    const std::vector<double> &boundary_condition,
    const int                  check) // whether to check the order of t
    :
    _t(t),
    _y(y), _method(method), _boundary_condition(boundary_condition) {
  poly = PPoly<double>();
  if (check) {
    if (!std::is_sorted(t.begin(), t.end())) {
      std::vector<int> idx(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        idx[i] = i;
      }
      std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });
      // std::vector<Real> t_sorted(t.size()), y_sorted(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        _t[i] = t[idx[i]];
        _y[i] = y[idx[i]];
      }
    }
  }
  interpolate(t, y, method, boundary_condition); // interpolate
}

template <>
void PPInterpolate<3, double>::interpolate(
    const std::vector<double> &t, // nodes
    const std::vector<double> &y, // values
    const int method,             // 0 for periodic, 1 for complete, 2 for natural, 3 for
                                  // not-a-knot.
    const std::vector<double> &boundary_condition) // boundary condition
{
  switch (method) {
  case 0: {
    int                 t_size = t.size();
    std::vector<double> K_i(t_size - 1, 0);

    // Compute K_i = f[x_i, x_{i+1}]

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    Eigen::SparseMatrix<double> A(t_size, t_size);
    Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);

    std::vector<Eigen::Triplet<double>>
        triplets; // For batch insertion into the sparse matrix
    triplets.reserve(4 * (t_size - 2));

    // Precompute mu_i and lambda_i to avoid redundant calculations

    for (int i = 0; i < t_size - 2; ++i) {
      double dt_i = t[i + 2] - t[i + 1];
      double dt_im1 = t[i + 1] - t[i];
      double mu_i = dt_im1 / dt_i;
      double lambda_i = dt_i / dt_im1;

      b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);

      // Store the triplets for A matrix in batch
      // Using thread-safe data structure (e.g., critical section or atomic
      // operations) may be required here

      {
        triplets.push_back(Eigen::Triplet<double>(i, i, lambda_i));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, 2.0));
        triplets.push_back(Eigen::Triplet<double>(i, i + 2, mu_i));
      }
    }

    // Boundary condition handling for the last row
    triplets.push_back(Eigen::Triplet<double>(t_size - 2, 0, 1.0));
    triplets.push_back(Eigen::Triplet<double>(t_size - 2, t_size - 1, -1.0));
    b(t_size - 2) = 0.0;

    // Boundary condition for the last element in A and b
    double dt_0 = t[1] - t[0];
    double dt_n = t[t_size - 1] - t[t_size - 2];
    triplets.push_back(Eigen::Triplet<double>(t_size - 1, 0, 4.0 / dt_0));
    triplets.push_back(Eigen::Triplet<double>(t_size - 1, 1, 2.0 / dt_0));
    triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 2, 2.0 / dt_n));
    triplets.push_back(Eigen::Triplet<double>(t_size - 1, t_size - 1, 4.0 / dt_n));

    b(t_size - 1) = (6.0 * K_i[t_size - 2] / dt_n + 6.0 * K_i[0] / dt_0);

    // Fill matrix A in batch
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXd c = solver.solve(b);

    // Prepare the cubic spline coefficients
    std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(4, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      double tmp = t[i + 1] - t[i];
      coeffs[i][0] = y[i];
      coeffs[i][1] = c(i);
      coeffs[i][2] = (3.0 * K_i[i] - 2.0 * c(i) - c(i + 1)) / tmp;
      coeffs[i][3] = (c(i) - 2.0 * K_i[i] + c(i + 1)) / (tmp * tmp);
    }
    poly = PPoly<double>(coeffs, t, 0); // Interpolate
    return;
  } break;

  case 1: {
    int                 t_size = t.size();
    std::vector<double> K_i(t_size - 1, 0);

    // compute K_i = f[x_i, x_{i+1}]

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    // init SparseMatrix A and Vector b
    Eigen::SparseMatrix<double>         A(t_size - 2, t_size - 2);
    Eigen::VectorXd                     b = Eigen::VectorXd::Zero(t_size - 2);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * (t_size - 2));

    // 填充矩阵 A 和向量 b

    for (int i = 0; i < t_size - 2; ++i) {
      double mu_i = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
      double lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
      b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);

      {
        triplets.emplace_back(i, i, 2.0);
        if (i == 0) {
          triplets.emplace_back(i, i + 1, mu_i);
        } else if (i < t_size - 3) {
          triplets.emplace_back(i, i - 1, mu_i);
          triplets.emplace_back(i, i + 1, lambda_i);
        } else {
          triplets.emplace_back(i, i - 1, lambda_i);
        }
      }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());

    // 使用稀疏矩阵求解器
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXd c = solver.solve(b);

    // 构建多项式系数
    std::vector<std::vector<double>> coeffs(t_size - 1, std::vector<double>(4, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      double tmp = t[i + 1] - t[i];
      coeffs[i][0] = y[i];

      if (i == 0) {
        coeffs[i][1] = boundary_condition[0];
        coeffs[i][2] = (3 * K_i[i] - 2 * boundary_condition[0] - c(0)) / tmp;
        coeffs[i][3] = (boundary_condition[0] - 2 * K_i[i] + c(0)) / (tmp * tmp);
      } else if (i < t_size - 2) {
        coeffs[i][1] = c(i - 1);
        coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - c(i)) / tmp;
        coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + c(i)) / (tmp * tmp);
      } else {
        coeffs[i][1] = c(i - 1);
        coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - boundary_condition[1]) / tmp;
        coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + boundary_condition[1]) / (tmp * tmp);
      }
    }

    // 创建 PPoly 对象
    poly = PPoly<double>(coeffs, t, 0);
    return;
  } break;

  case 2: {
    /**
     * Natural Cubic Spline
     */
    int                 t_size = t.size();
    std::vector<double> K_i(t_size - 1, 0);

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }
    std::vector<double> J_i(t_size - 2, 0);

    for (int i = 0; i < t_size - 2; ++i) {
      J_i[i] = (K_i[i + 1] - K_i[i]) / (t[i + 2] - t[i]);
    }

    Eigen::VectorXd b = Eigen::VectorXd::Zero(t_size - 2);
    // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(t_size - 2, t_size - 2);
    Eigen::SparseMatrix<double>         A(t_size - 2, t_size - 2);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * (t_size - 2));

    for (int i = 0; i < t_size - 2; ++i) {
      b(i) = 6 * J_i[i];

      double mu_i = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
      double lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
      // A(i,i) = 2;

      {
        triplets.emplace_back(i, i, 2.0);
        if (i == 0) {
          // A(i, i+1) = lambda_i;
          triplets.emplace_back(i, i + 1, lambda_i);
        } else if (i < t_size - 3) {
          // A(i, i-1) = mu_i;
          // A(i, i+1) = lambda_i;
          triplets.emplace_back(i, i - 1, mu_i);
          triplets.emplace_back(i, i + 1, lambda_i);
        } else {
          // A(i, i-1) = mu_i;
          triplets.emplace_back(i, i - 1, mu_i);
        }
      }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    // solve Ax = b

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXd c = solver.solve(b);

    std::vector<std::vector<double>> coeffs(t.size() - 1, std::vector<double>(4, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      coeffs[i][0] = y[i];
      if (i == 0) {

        double tmp = t[i + 1] - t[i];
        coeffs[i][1] = K_i[i] - 1.0 / 6.0 * c(0) * tmp;
        coeffs[i][2] = 0;
        coeffs[i][3] = c(0) / (6.0 * tmp);
      } else if (i < t_size - 2) {
        double tmp = t[i + 1] - t[i];
        coeffs[i][1] = K_i[i] - 1.0 / 6.0 * (c(i) + 2 * c(i - 1)) * tmp;
        coeffs[i][2] = c(i - 1) / 2.0;
        coeffs[i][3] = (c(i) - c(i - 1)) / (6.0 * tmp);
      } else {
        double tmp = t[i + 1] - t[i];
        coeffs[i][1] = K_i[i] - 1.0 / 3.0 * c(i - 1) * tmp;
        coeffs[i][2] = c(i - 1) / 2.0;
        coeffs[i][3] = -c(i - 1) / (6.0 * tmp);
      }
    }
    poly = PPoly<double>(coeffs, t, 0); // interpolate
    return;
  } break;

  default:
    throw std::invalid_argument("method must be 0,1,2");
    return;
    break;
  }
  return;
}

template <>
PPInterpolate<1, mpf_class>::PPInterpolate(
    const std::vector<mpf_class> &t, // nodes
    const std::vector<mpf_class> &y, // values
    const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                      // not-a-knot.
    const std::vector<mpf_class> &boundary_condition,
    const int                     check) // whether to check the order of t
    :
    _t(t),
    _y(y), _method(method), _boundary_condition(boundary_condition) {
  poly = PPoly<mpf_class>();
  if (check) {
    if (!std::is_sorted(t.begin(), t.end())) {
      std::vector<int> idx(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        idx[i] = i;
      }
      std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });
      // std::vector<Real> t_sorted(t.size()), y_sorted(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        _t[i] = t[idx[i]];
        _y[i] = y[idx[i]];
      }
    }
  }
  interpolate(t, y, method, boundary_condition); // interpolate
}

template <>
void PPInterpolate<1, mpf_class>::interpolate(
    const std::vector<mpf_class> &t, // nodes
    const std::vector<mpf_class> &y, // values
    const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                      // not-a-knot.
    const std::vector<mpf_class> &boundary_condition) // boundary condition
{
  // unused variable, how to solve it.
  (void)method;
  (void)boundary_condition;
  std::vector<std::vector<mpf_class>> A(t.size() - 1, std::vector<mpf_class>(2, 0));

  for (size_t i = 0; i < t.size() - 1; ++i) {
    A[i][0] = y[i];
    A[i][1] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
  }
  poly = PPoly<mpf_class>(A, t, 0); // interpolate
}

template <>
PPInterpolate<2, mpf_class>::PPInterpolate(
    const std::vector<mpf_class> &t, // nodes
    const std::vector<mpf_class> &y, // values
    const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                      // not-a-knot.
    const std::vector<mpf_class> &boundary_condition,
    const int                     check) // whether to check the order of t
    :
    _t(t),
    _y(y), _method(method), _boundary_condition(boundary_condition) {
  poly = PPoly<mpf_class>();
  if (check) {
    if (!std::is_sorted(t.begin(), t.end())) {
      std::vector<int> idx(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        idx[i] = i;
      }
      std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });
      // std::vector<Real> t_sorted(t.size()), y_sorted(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        _t[i] = t[idx[i]];
        _y[i] = y[idx[i]];
      }
    }
  }
  interpolate(t, y, method, boundary_condition); // interpolate
}

template <>
void PPInterpolate<2, mpf_class>::interpolate(
    const std::vector<mpf_class> &t, // nodes
    const std::vector<mpf_class> &y, // values
    const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                      // not-a-knot.
    const std::vector<mpf_class> &boundary_condition) // boundary condition
{
  switch (method) {
  case 0: {
    int                    t_size = t.size();
    std::vector<mpf_class> K_i(t_size - 1, 0);

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    Eigen::SparseMatrix<mpf_class> A(t_size, t_size);
    // Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size);

    std::vector<Eigen::Triplet<mpf_class>>
        triplets; // For batch insertion into the sparse matrix

    for (int i = 0; i < t_size - 1; i++) {
      b(i) = 2 * K_i[i];

      {
        triplets.push_back(Eigen::Triplet<mpf_class>(i, i, 1));
        triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, 1));
      }
    }
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, 0, 1));
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 1, -1));
    b(t_size - 1) = 0;
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    // Eigen::VectorXd c = solver.solve(b);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

    std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(3, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      coeffs[i][0] = y[i];
      coeffs[i][1] = c(i);
      coeffs[i][2] = (c(i + 1) - c(i)) / (2 * (t[i + 1] - t[i]));
    }
    poly = PPoly<mpf_class>(coeffs, t, 0); // interpolate
    return;
  } break;

  case 1: {
    int                    t_size = t.size();
    std::vector<mpf_class> K_i(t_size - 1, 0);

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    std::vector<mpf_class> c(t_size, 0);
    c[0] = boundary_condition[0];

    for (int i = 1; i < t_size; ++i) {
      c[i] = 2 * K_i[i - 1] - c[i - 1];
    }

    std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(3, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      coeffs[i][0] = y[i];
      coeffs[i][1] = c[i];
      coeffs[i][2] = (c[i + 1] - c[i]) / (2 * (t[i + 1] - t[i]));
    }
    poly = PPoly<mpf_class>(coeffs, t, 0); // interpolate
    return;
  } break;

  default:
    throw std::invalid_argument("method must be 0,1");
    return;
    break;
  }
  return;
}

template <>
PPInterpolate<3, mpf_class>::PPInterpolate(
    const std::vector<mpf_class> &t, // nodes
    const std::vector<mpf_class> &y, // values
    const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                      // not-a-knot.
    const std::vector<mpf_class> &boundary_condition,
    const int                     check) // whether to check the order of t
    :
    _t(t),
    _y(y), _method(method), _boundary_condition(boundary_condition) {
  poly = PPoly<mpf_class>();
  if (check) {
    if (!std::is_sorted(t.begin(), t.end())) {
      std::vector<int> idx(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        idx[i] = i;
      }
      std::sort(idx.begin(), idx.end(), [&](int i, int j) { return t[i] < t[j]; });
      // std::vector<Real> t_sorted(t.size()), y_sorted(t.size());

      for (size_t i = 0; i < t.size(); ++i) {
        _t[i] = t[idx[i]];
        _y[i] = y[idx[i]];
      }
    }
  }
  interpolate(t, y, method, boundary_condition); // interpolate
}

template <>
void PPInterpolate<3, mpf_class>::interpolate(
    const std::vector<mpf_class> &t, // nodes
    const std::vector<mpf_class> &y, // values
    const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                      // not-a-knot.
    const std::vector<mpf_class> &boundary_condition) // boundary condition
{
  switch (method) {
  case 0: {
    int                    t_size = t.size();
    std::vector<mpf_class> K_i(t_size - 1, 0);

    // Compute K_i = f[x_i, x_{i+1}]

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    Eigen::SparseMatrix<mpf_class> A(t_size, t_size);
    // Eigen::VectorXd             b = Eigen::VectorXd::Zero(t_size);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size);

    std::vector<Eigen::Triplet<mpf_class>>
        triplets; // For batch insertion into the sparse matrix
    triplets.reserve(4 * (t_size - 2));

    // Precompute mu_i and lambda_i to avoid redundant calculations

    for (int i = 0; i < t_size - 2; ++i) {
      mpf_class dt_i = t[i + 2] - t[i + 1];
      mpf_class dt_im1 = t[i + 1] - t[i];
      mpf_class mu_i = dt_im1 / dt_i;
      mpf_class lambda_i = dt_i / dt_im1;

      b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);

      // Store the triplets for A matrix in batch
      // Using thread-safe data structure (e.g., critical section or atomic
      // operations) may be required here

      {
        triplets.push_back(Eigen::Triplet<mpf_class>(i, i, lambda_i));
        triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 1, 2.0));
        triplets.push_back(Eigen::Triplet<mpf_class>(i, i + 2, mu_i));
      }
    }

    // Boundary condition handling for the last row
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 2, 0, 1.0));
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 2, t_size - 1, -1.0));
    b(t_size - 2) = 0.0;

    // Boundary condition for the last element in A and b
    mpf_class dt_0 = t[1] - t[0];
    mpf_class dt_n = t[t_size - 1] - t[t_size - 2];
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, 0, 4.0 / dt_0));
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, 1, 2.0 / dt_0));
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 2, 2.0 / dt_n));
    triplets.push_back(Eigen::Triplet<mpf_class>(t_size - 1, t_size - 1, 4.0 / dt_n));

    b(t_size - 1) = (6.0 * K_i[t_size - 2] / dt_n + 6.0 * K_i[0] / dt_0);

    // Fill matrix A in batch
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    // Eigen::VectorXd c = solver.solve(b);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

    // Prepare the cubic spline coefficients
    std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(4, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      mpf_class tmp = t[i + 1] - t[i];
      coeffs[i][0] = y[i];
      coeffs[i][1] = c(i);
      coeffs[i][2] = (3.0 * K_i[i] - 2.0 * c(i) - c(i + 1)) / tmp;
      coeffs[i][3] = (c(i) - 2.0 * K_i[i] + c(i + 1)) / (tmp * tmp);
    }
    poly = PPoly<mpf_class>(coeffs, t, 0); // Interpolate
    return;
  } break;

  case 1: {
    int                    t_size = t.size();
    std::vector<mpf_class> K_i(t_size - 1, 0);

    // compute K_i = f[x_i, x_{i+1}]

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }

    // init SparseMatrix A and Vector b
    Eigen::SparseMatrix<mpf_class> A(t_size - 2, t_size - 2);
    // Eigen::VectorXd                     b = Eigen::VectorXd::Zero(t_size - 2);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size - 2);
    std::vector<Eigen::Triplet<mpf_class>>      triplets;
    triplets.reserve(3 * (t_size - 2));

    // 填充矩阵 A 和向量 b

    for (int i = 0; i < t_size - 2; ++i) {
      mpf_class mu_i = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
      mpf_class lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
      b(i) = 3 * (mu_i * K_i[i + 1] + lambda_i * K_i[i]);

      {
        triplets.emplace_back(i, i, 2.0);
        if (i == 0) {
          triplets.emplace_back(i, i + 1, mu_i);
        } else if (i < t_size - 3) {
          triplets.emplace_back(i, i - 1, mu_i);
          triplets.emplace_back(i, i + 1, lambda_i);
        } else {
          triplets.emplace_back(i, i - 1, lambda_i);
        }
      }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());

    // 使用稀疏矩阵求解器
    Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    // Eigen::VectorXd c = solver.solve(b);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

    // 构建多项式系数
    std::vector<std::vector<mpf_class>> coeffs(t_size - 1, std::vector<mpf_class>(4, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      mpf_class tmp = t[i + 1] - t[i];
      coeffs[i][0] = y[i];

      if (i == 0) {
        coeffs[i][1] = boundary_condition[0];
        coeffs[i][2] = (3 * K_i[i] - 2 * boundary_condition[0] - c(0)) / tmp;
        coeffs[i][3] = (boundary_condition[0] - 2 * K_i[i] + c(0)) / (tmp * tmp);
      } else if (i < t_size - 2) {
        coeffs[i][1] = c(i - 1);
        coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - c(i)) / tmp;
        coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + c(i)) / (tmp * tmp);
      } else {
        coeffs[i][1] = c(i - 1);
        coeffs[i][2] = (3 * K_i[i] - 2 * c(i - 1) - boundary_condition[1]) / tmp;
        coeffs[i][3] = (c(i - 1) - 2 * K_i[i] + boundary_condition[1]) / (tmp * tmp);
      }
    }

    // 创建 PPoly 对象
    poly = PPoly<mpf_class>(coeffs, t, 0);
    return;
  } break;

  case 2: {
    /**
     * Natural Cubic Spline
     */
    int                    t_size = t.size();
    std::vector<mpf_class> K_i(t_size - 1, 0);

    for (int i = 0; i < t_size - 1; ++i) {
      K_i[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]);
    }
    std::vector<mpf_class> J_i(t_size - 2, 0);

    for (int i = 0; i < t_size - 2; ++i) {
      J_i[i] = (K_i[i + 1] - K_i[i]) / (t[i + 2] - t[i]);
    }

    // Eigen::VectorXd b = Eigen::VectorXd::Zero(t_size - 2);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> b(t_size - 2);
    // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(t_size - 2, t_size - 2);
    Eigen::SparseMatrix<mpf_class>         A(t_size - 2, t_size - 2);
    std::vector<Eigen::Triplet<mpf_class>> triplets;
    triplets.reserve(3 * (t_size - 2));

    for (int i = 0; i < t_size - 2; ++i) {
      b(i) = 6 * J_i[i];

      mpf_class mu_i = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
      mpf_class lambda_i = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
      // A(i,i) = 2;

      {
        triplets.emplace_back(i, i, 2.0);
        if (i == 0) {
          // A(i, i+1) = lambda_i;
          triplets.emplace_back(i, i + 1, lambda_i);
        } else if (i < t_size - 3) {
          // A(i, i-1) = mu_i;
          // A(i, i+1) = lambda_i;
          triplets.emplace_back(i, i - 1, mu_i);
          triplets.emplace_back(i, i + 1, lambda_i);
        } else {
          // A(i, i-1) = mu_i;
          triplets.emplace_back(i, i - 1, mu_i);
        }
      }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    // solve Ax = b

    Eigen::SparseLU<Eigen::SparseMatrix<mpf_class>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    // Eigen::VectorXd c = solver.solve(b);
    Eigen::Matrix<mpf_class, Eigen::Dynamic, 1> c = solver.solve(b);

    std::vector<std::vector<mpf_class>> coeffs(
        t.size() - 1, std::vector<mpf_class>(4, 0));

    for (int i = 0; i < t_size - 1; ++i) {
      coeffs[i][0] = y[i];
      if (i == 0) {

        mpf_class tmp = t[i + 1] - t[i];
        coeffs[i][1] = K_i[i] - 1.0 / 6.0 * c(0) * tmp;
        coeffs[i][2] = 0;
        coeffs[i][3] = c(0) / (6.0 * tmp);
      } else if (i < t_size - 2) {
        mpf_class tmp = t[i + 1] - t[i];
        coeffs[i][1] = K_i[i] - 1.0 / 6.0 * (c(i) + 2 * c(i - 1)) * tmp;
        coeffs[i][2] = c(i - 1) / 2.0;
        coeffs[i][3] = (c(i) - c(i - 1)) / (6.0 * tmp);
      } else {
        mpf_class tmp = t[i + 1] - t[i];
        coeffs[i][1] = K_i[i] - 1.0 / 3.0 * c(i - 1) * tmp;
        coeffs[i][2] = c(i - 1) / 2.0;
        coeffs[i][3] = -c(i - 1) / (6.0 * tmp);
      }
    }
    poly = PPoly<mpf_class>(coeffs, t, 0); // interpolate
    return;
  } break;

  default:
    throw std::invalid_argument("method must be 0,1,2");
    return;
    break;
  }
  return;
}

template class PPInterpolate<1, double>;
template class PPInterpolate<2, double>;
template class PPInterpolate<3, double>;
template class PPInterpolate<1, mpf_class>;
template class PPInterpolate<2, mpf_class>;
template class PPInterpolate<3, mpf_class>;