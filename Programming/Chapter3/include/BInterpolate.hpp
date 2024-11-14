/**
 * @file BInterpolate.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Use B-spline to interpolate the data
 * @version 0.1
 * @date 2024-11-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */


#ifndef BINTERPOLATE_HPP
#define BINTERPOLATE_HPP

#include "BSpline.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include <vector>
#include <cmath>
#include <string>
#include <omp.h>

template <int N>
class BInterpolate
{
private:
    /**
     * nodes of the interpolation
     */
    std::vector<double> _t;

    /**
     * interpolation data
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
    BSpline _bspline;

public:
    BInterpolate() = default;

    BInterpolate(const std::vector<double> &t, 
                 const std::vector<double> &y, 
                 const int &method=0, 
                 const std::vector<double> &boundary_condition = std::vector<double>(N, 0.0),
                 const int check=0);

    void
    interpolate(const std::vector<double> &t, 
                 const std::vector<double> &y, 
                 const int &method, 
                 const std::vector<double> &boundary_condition = std::vector<double>(N, 0.0));
    
    BSpline
    getBSpline() const;

    double
    operator()(const double x);

    double 
    derivative(const double x, const int n);
};




#endif // BINTERPOLATE_HPP