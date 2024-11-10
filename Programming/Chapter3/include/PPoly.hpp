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

#include <iostream>
#include <vector>
#include <cmath>

/**
 * @brief The definition of the class PPoly, which is a polynomial class.
 * 
 * @tparam N The order of the polynomial, \f$ \mathbb{S}_{N}^{N-1} \f$
 * 
 * @details The class PPoly is a polynomial class, which is used to represent a piecewise polynomial function.
 * In each interval \f$ [t_{i}, t_{i+1}] \f$, the polynomial is represented as
 * 
 *  
 */
class PPoly 
{
private:
    /**
     * The coefficients of pp-form
     */
    std::vector<std::vector<double>> _coeffs;

    /**
     * The knots of pp-form
     */
    std::vector<double> _t;

public:
    /**
     * Construct a new PPoly object
     */
    PPoly();
    

    /**
     * The function use the coefficients and knots to construct a new PPoly object.
     */
    PPoly(const std::vector<std::vector<double>> &coeffs, /** The coefficients of pp-form */ 
          const std::vector<double> &t, /** The knots of pp-form */
          const int check = 1); 
    
    /**
     * @brief The function return the value of the polynomial at a given point.
     * 
     * @param x the x value
     * @return double \f$ p(x) \f$
     */
    double 
    operator()(double x) const;
    

    /**
     * @brief The function return the value of the derivative of the polynomial at a given point.
     * 
     * @param x the x value
     * @param n the order of derivative
     * @return double \f$ p^{(n)}(x) \f$
     */
    double 
    derivative(double x, int n = 1) const;
    
    /**
     * @brief The function return the value of the integral of the polynomial at a given interval.
     * 
     * @param a The left boundary of the integral
     * @param b The right boundary of the integral
     * @return double \f$ \int_{a}^{b} p(x) dx \f$
     */
    double 
    integral(double a, double b) const;
    
    // /**
    //  * @brief The function use the coefficients and knots to construct a new PPoly object.
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
     * @brief The function return the index of the interval that contains the given point.
     * 
     * @param x The point
     * @return int 
     */
    int 
    findInterval(double x) const;

    // overload =
    PPoly &
    operator=(const PPoly &other);

};



#endif // PPOLY_HPP