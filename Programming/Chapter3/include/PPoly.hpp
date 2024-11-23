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
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <vector>

#define EIGEN_USE_THREADS

/**
 * @brief The definition of the class PPoly, which is a polynomial class.
 *
 * @tparam N The order of the polynomial, \f$ \mathbb{S}_{N}^{N-1} \f$
 *
 * @details The class PPoly is a polynomial class, which is used to represent a
 * piecewise polynomial function. In each interval \f$ [t_{i}, t_{i+1}] \f$, the
 * polynomial is represented as
 *
 *
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
     * Construct a new PPoly object
     */
    PPoly();

    /**
     * The function use the coefficients and knots to construct a new PPoly
     * object.
     */
    PPoly(const std::vector<std::vector<Real>>
                                  &coeffs, /** The coefficients of pp-form */
          const std::vector<Real> &t,      /** The knots of pp-form */
          const int                check = 0);

    /**
     * @brief Construct a new PPoly object
     *
     * @param poly
     */
    PPoly(const PPoly &poly);

    /**
     * @brief The function return the value of the polynomial at a given point.
     *
     * @param x the x value
     * @return Real \f$ p(x) \f$
     */
    Real
    operator()(Real x) const;

    /**
     * @brief The function return the value of the derivative of the polynomial
     * at a given point.
     *
     * @param x the x value
     * @param n the order of derivative
     * @return Real \f$ p^{(n)}(x) \f$
     */
    Real
    derivative(Real x, int n = 1) const;

    /**
     * @brief The function return the value of the integral of the polynomial at
     * a given interval.
     *
     * @param a The left boundary of the integral
     * @param b The right boundary of the integral
     * @return Real \f$ \int_{a}^{b} p(x) dx \f$
     */
    Real
    integral(Real a, Real b) const;

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
     * @brief The function return the index of the interval that contains the
     * given point.
     *
     * @param x The point
     * @return int
     */
    int
    findInterval(Real x) const;

    // overload =
    PPoly &
    operator=(const PPoly &other);
};

#include "PPoly.tpp"

#endif // PPOLY_HPP