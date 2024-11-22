#include "BInterpolate.hpp"

template <>
void BInterpolate<1, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int                 &method,
    const std::vector<double> &boundary_condition) {
  (void)method;
  (void)boundary_condition;
  _bspline = BSpline(y, t, 1); // default method is 1
  // std::cout << "Interpolation finished." << std::endl;
}

template <>
void BInterpolate<3, double>::interpolate(
    const std::vector<double> &t,
    const std::vector<double> &y,
    const int &method, // 0 for periodic, 1 for complete, 2 for natural, 3 for
                       // not-a-knot.
    const std::vector<double> &boundary_condition) {
  if (method == 0) {
    int                         t_size = t.size();
    Eigen::SparseMatrix<double> A(t_size + 2, t_size + 2);
    // Eigen::VectorXd b(t_size+2);
    Eigen::Matrix<double, Eigen::Dynamic, 1> b(t_size + 2);
    std::vector<Eigen::Triplet<double>>      triplets;
    triplets.reserve(4 * (t_size + 2));
    std::vector<double> tmp_coeff(t_size + 2, 1.0);
    BSpline             tmp_spline(tmp_coeff, t, 3);
    for (int i = 0; i < t_size; ++i) {
      std::vector<double> basis = tmp_spline.get_basis(t[i]);
      if (i == 0) {
        triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[2]));
      } else {
        triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[3]));
      }
      b(i) = y[i];
    }
    std::vector<double> diff_basis1 = tmp_spline.basis_derivative(t[0], 1);
    std::vector<double> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 1);
    triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
    triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis1[2]));
    triplets.push_back(Eigen::Triplet<double>(t_size, t_size - 1, -diff_basis2[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size, t_size, -diff_basis2[2]));
    triplets.push_back(Eigen::Triplet<double>(t_size, t_size + 1, -diff_basis2[3]));
    b(t_size) = 0.0;

    diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
    diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, 0, diff_basis1[0]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, 1, diff_basis1[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, 2, diff_basis1[2]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size - 1, -diff_basis2[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size, -diff_basis2[2]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size + 1, -diff_basis2[3]));
    b(t_size + 1) = 0.0;

    A.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    // Eigen::VectorXd x = solver.solve(b);
    Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(b);

    std::vector<double> coeff(x.data(), x.data() + x.size());
    _bspline = BSpline(coeff, t, 3);
    return;
  } else if (method == 1) {
    int                         t_size = t.size();
    Eigen::SparseMatrix<double> A(t_size + 2, t_size + 2);
    // Eigen::VectorXd b(t_size+2);
    Eigen::Matrix<double, Eigen::Dynamic, 1> b(t_size + 2);
    std::vector<Eigen::Triplet<double>>      triplets;
    triplets.reserve(4 * (t_size + 2));
    std::vector<double> tmp_coeff(t_size + 2, 1.0);
    BSpline             tmp_spline(tmp_coeff, t, 3);
    for (int i = 0; i < t_size; ++i) {
      std::vector<double> basis = tmp_spline.get_basis(t[i]);
      if (i == 0) {
        triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[2]));
      } else {
        triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[3]));
      }
      b(i) = y[i];
    }

    std::vector<double> diff_basis = tmp_spline.basis_derivative(t[0], 1);
    triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis[0]));
    triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis[2]));
    b(t_size) = boundary_condition[0];

    diff_basis = tmp_spline.basis_derivative(t[t_size - 1], 1);
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size - 1, diff_basis[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size, diff_basis[2]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size + 1, diff_basis[3]));
    b(t_size + 1) = boundary_condition[1];
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    // Eigen::VectorXd x = solver.solve(b);
    Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(b);
    std::vector<double>                      coeffs(x.data(), x.data() + x.size());

    // // convert x to std::vector
    // tmp_coeff = std::vector<double>(x.data(), x.data() + x.size());

    // _bspline = BSpline(tmp_coeff, t, 3);
    _bspline = BSpline(coeffs, t, 3);
    return;

  } else if (method == 2) {
    int                                      t_size = t.size();
    Eigen::SparseMatrix<double>              A(t_size + 2, t_size + 2);
    Eigen::Matrix<double, Eigen::Dynamic, 1> b(t_size + 2);
    std::vector<Eigen::Triplet<double>>      triplets;
    triplets.reserve(4 * (t_size + 2));
    std::vector<double> tmp_coeff(t_size + 2, 1.0);
    BSpline             tmp_spline(tmp_coeff, t, 3);
    for (int i = 0; i < t_size; ++i) {
      std::vector<double> basis = tmp_spline.get_basis(t[i]);
      if (i == 0) {
        triplets.push_back(Eigen::Triplet<double>(i, i, basis[0]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[1]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[2]));
      } else {
        triplets.push_back(Eigen::Triplet<double>(i, i, basis[1]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, basis[2]));
        triplets.push_back(Eigen::Triplet<double>(i, i + 2, basis[3]));
      }
      b(i) = y[i];
    }

    std::vector<double> diff_basis1 = tmp_spline.basis_derivative(t[0], 2);
    std::vector<double> diff_basis2 = tmp_spline.basis_derivative(t[t_size - 1], 2);
    triplets.push_back(Eigen::Triplet<double>(t_size, 0, diff_basis1[0]));
    triplets.push_back(Eigen::Triplet<double>(t_size, 1, diff_basis1[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size, 2, diff_basis1[2]));
    b(t_size) = 0.0;
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size - 1, diff_basis2[1]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size, diff_basis2[2]));
    triplets.push_back(Eigen::Triplet<double>(t_size + 1, t_size + 1, diff_basis2[3]));
    b(t_size + 1) = 0.0;

    A.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    // Eigen::VectorXd x = solver.solve(b);
    Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(b);

    std::vector<double> coeff(x.data(), x.data() + x.size());
    _bspline = BSpline(coeff, t, 3);
    return;
  } else if (method == 4) {
    // to be done
  } else {
    throw std::invalid_argument("method must be 0, 1, 2");
  }
}

template class BInterpolate<1, double>;
template class BInterpolate<3, double>;