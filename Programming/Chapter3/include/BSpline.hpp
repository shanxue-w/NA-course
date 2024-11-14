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

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <omp.h>
#include <lapacke.h>
#include <lapack.h>

#define EIGEN_USE_LAPACKE
#define EIGEN_USE_LAPACKE_STRICT

class BSpline
{
    /**
     * @brief The Bspline is in the form of 
     * \f[ B(x) = \Sum_{j=i-n+1}^{i+1} c_j B_j^n(x), ~ x\in[x_i, x_{i+1}] \f]
     * So the computation of B_j^n(x) is very important.
     */
private:
    std::vector<double> _coeffs;
    std::vector<double> _t;
    std::vector<double> _basis;
    std::vector<std::vector<double>> _total_basis;
    std::vector<double> _derivative_basis;
    int _n;

public:
    BSpline() = default;

    BSpline(const std::vector<double> coeffs, 
            const std::vector<double> t, 
            const int                 n,
            const int                 check=0);
    
    BSpline(const BSpline& other);

    int
    get_interval(const double x) const;

    BSpline& 
    operator=(const BSpline& other);

    double 
    operator()(const double x);

    std::vector<double> 
    get_coeffs() const;

    std::vector<double> 
    get_t() const;

    int 
    get_n() const;

    std::vector<double>
    get_basis(const double x);

    void
    get_total_basis(const double x);

    double
    cal_derivative_basis(const int interval,
                         const int j,
                         const int degree,
                         const int n);

    std::vector<double>
    basis_derivative(const double x, 
                     const int    n);

    double
    derivative(const double x, 
               const int    n);

};





#endif // BSPLINE_HPP