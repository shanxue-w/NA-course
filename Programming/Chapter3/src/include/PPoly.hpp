/**
 * @file PPoly.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The definition of the class PPoly, which is a polynomial class.
 * @version 0.1
 * @date 2024-11-09
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef PPOLY_HPP
#define PPOLY_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <omp.h>
#include <vector>

#define EIGEN_USE_THREADS
#define EIGEN_USE_LAPACKE

/**
 * @brief The definition of the class PPoly, which is a polynomial class.
 *
 * @tparam N The order of the polynomial, \f$ \mathbb{S}_{N}^{N-1} \f$
 *
 * @details The class PPoly is a polynomial class, which is used to represent a
 * piecewise polynomial function. In each interval \f$ [t_{i}, t_{i+1}] \f$, the
 * polynomial is represented as
 * \f[ p_{i}(x) = \sum_{j=0}^{N-1} \alpha_{ij} (x - t_{i})^{j}, x \in [t_{i}, t_{i+1}] \f].
 *
 * The coefficients \f$ \alpha_{ij} \f$ are stored in a matrix, and the knots \f$ t_{i} \f$ are stored in a vector.
 * The class provides a constructor to initialize the polynomial, and an operator() to evaluate the polynomial at a given point.
 * The class also provides a method to compute the derivative of the polynomial.
 */
template <typename Real = double> // the type of the coefficients
class PPoly
{
private:
    /**
     * The coefficients of pp-form
     */
    std::vector<std::vector<Real>> _coeffs;
    // Eigen::MatrixXd _coeffs;

    /**
     * The knots of pp-form
     */
    std::vector<Real> _t;

public:
    /**
     * @brief Construct a new PPoly<Real>::PPoly object
     *
     * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
     *
     * @details The default constructor of the class PPoly.
     * The coefficients and knots are initialized as empty.
     */
    PPoly();

    /**
     * @brief Construct a new PPoly<Real>::PPoly object
     *
     * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
     * @param coeffs The coefficients of pp-form
     * @param t The knots of pp-form
     * @param check Whether to check the input data, default is 0.
     *
     * @details The pp-form polynomial has the form
     * \f[
     * p_i(x) = \sum_{j=0}^{n} a_{i,j} (x - t_i)^j, x \in [t_i, t_{i+1}].
     * \f]
     * \f$a_{i,j}\f$ are the `coeffs`, and \f$t_i\f$ are the `t`.
     *
     *
     * make sure the t is sorted, and the size of t is equal to the size of coeffs.
     * coeffs is corresponding to the interval [t[i], t[i+1]],
     * so when sorting t, the coeffs should be sorted too.
     */
    PPoly(const std::vector<std::vector<Real>> &coeffs, /** The coefficients of pp-form */
          const std::vector<Real>              &t,      /** The knots of pp-form */
          const int                             check = 0);

    /**
     * @brief Construct a new PPoly<Real>::PPoly object
     *
     * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
     * @param other The object to copy, class `PPoly`.
     */
    PPoly(const PPoly &poly);

    /**
     * @brief Overload the operator() to calculate the value of the polynomial at x.
     *
     * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
     * @param x The value to calculate the polynomial.
     * @return Real The value of the polynomial at x.
     *
     * @details The function return the value of the polynomial at a given point.
     *
     * First find the interval that x belongs to.
     * Then use the Horner's method to calculate the value of the polynomial.
     *
     * \f[
     * p_i(x) = \sum_{j=0}^{n} a_{i,j} (x - t_i)^j
     * = a_{i,0} + (x-t_i)(a_{i,1} + (x-t_i)(a_{i,2} + \cdots ))
     * \f]
     */
    Real
    operator()(Real x) const;

    /**
     * @brief Calculate the n-th derivative of the polynomial at x.
     *
     * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
     * @param x The value to calculate the derivative.
     * @param n The order of the derivative.
     * @return Real The value of the n-th derivative of the polynomial at x.
     *
     * @details If the derivative order n is greater than the degree of the polynomial,
     * return 0. If n is negative, return 0. If x is out of bounds, return 0.
     *
     * Then just derive the polynomial and use Horner's method to evaluate the n-th
     * derivative at x.
     */
    Real
    derivative(Real x, int n = 1) const;

    // /**
    //  * @brief The function return the value of the integral of the polynomial at
    //  * a given interval.
    //  *
    //  * @param a The left boundary of the integral
    //  * @param b The right boundary of the integral
    //  * @return Real \f$ \int_{a}^{b} p(x) dx \f$
    //  */
    // Real
    // integral(Real a, Real b) const;

    // /**
    //  * @brief The function use the coefficients and knots to construct a new
    //  PPoly object.
    //  * In the form of
    //  * \f{equation}{
    //  * p(x) = \sum_{i=0}^{N} c_{i} (x-t_{j})^{i}
    //  * \f}
    //  *
    //  * @param os The output stream
    //  * @return std::ostream&
    //  */
    // std::ostream &operator<<(std::ostream &os) const;

    /**
     * @brief Find the interval that x belongs to.
     *
     * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
     * @param x The value to find the interval.
     * @return int
     *
     * @details When we input a value \f$ x \f$, we need to find the interval that
     * \f$ x \in[t_i, t_{i+1}) \f$, and return the index \f$ i \f$.
     * For the end point \f$ t_n \f$, we return \f$ n-1 \f$.
     *
     * If x is out of the range of the knots, return -1, means
     * we don't have the interval that x belongs to.
     *
     * @note For the knots already sorted, we can use binary search to find the
     * interval.
     *
     */
    int
    findInterval(Real x) const;

    /**
     * @brief Overload the operator= to copy the object.
     *
     * @tparam Real The type of data, can be `double`, `float`, `mpf_class` etc.
     * @param other The object to copy.
     * @return PPoly<Real>&
     */
    PPoly &
    operator=(const PPoly &other);
};

#include "PPoly.tpp"

#endif // PPOLY_HPP