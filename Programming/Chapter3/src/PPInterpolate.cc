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
#include <algorithm>
#include <iostream>

//
template <int N>
PPInterpolate<N>::PPInterpolate(
                  const std::vector<double> &t, // nodes
                  const std::vector<double> &y, // values
                  const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
                  const std::vector<double> &boundary_condition)
    : _method(method)
    , _boundary_condition(boundary_condition)
{
    poly = PPoly();
    _t = t;
    _y = y;
    if (!std::is_sorted(t.begin(), t.end()))
    {
        std::vector<int> idx(t.size());
        for (size_t i=0; i<t.size(); ++i)
            idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](int i, int j){return t[i] < t[j];});
        std::vector<double> t_sorted(t.size()), y_sorted(t.size());
        for (size_t i=0; i<t.size(); ++i)
        {
            t_sorted[i] = t[idx[i]];
            y_sorted[i] = y[idx[i]];
        }
        _t = t_sorted;
        _y = y_sorted;
    }
    interpolate(t, y, method, boundary_condition); // interpolate
}

// template <int N>
// void 
// PPInterpolate<N>::interpolate(
//                   const std::vector<double> &t, // nodes
//                   const std::vector<double> &y, // values
//                   const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
//                   const std::vector<double> &boundary_condition) // boundary condition
// {
//     (void) t;
//     (void) y;
//     (void) boundary_condition; // to be done
//     if (method == 0)
//     {
//         // To be done
//     }
//     else if (method == 1)
//     {
//         Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);

//         for (int i = 0; i < N; ++i) {
//             A(i, i) = 1.0;
//             A(0, i) = i + 1.0;
//         }

//         for (int i = 1; i < N; ++i) {
//             for (int j = i + 1; j < N; ++j) {
//                 A(i, j) = A(i, j - 1) + A(i - 1, j - 1);
//             }
//         }

        


//     }
// }

template <int N>
double
PPInterpolate<N>::operator()(double x) const
{
    return poly(x);
}



// template <>
// void PPInterpolate<1>::interpolate(
//                   const std::vector<double> &t, // nodes
//                   const std::vector<double> &y, // values
//                   const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
//                   const std::vector<double> &boundary_condition) // boundary condition
// {
//     // unused variable, how to solve it.
//     (void)method;
//     (void)boundary_condition;
//     std::vector<std::vector<double>> A(t.size()-1, std::vector<double>(2,0));
//     for (size_t i=0; i<t.size()-1; ++i)
//     {
//         A[i][0] = y[i];
//         A[i][1] = (y[i+1]-y[i])/(t[i+1]-t[i]);
//     }
//     poly = PPoly(A, t); // interpolate
// }

// template <>
// void PPInterpolate<2>::interpolate(
//                   const std::vector<double> &t, // nodes
//                   const std::vector<double> &y, // values
//                   const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
//                   const std::vector<double> &boundary_condition) // boundary condition
// {
//     (void)t;
//     (void)y;
//     (void)method;
//     (void)boundary_condition;
//     // To be done
//     return;
// }


// template <>
// PPInterpolate<3>::PPInterpolate(
//                   const std::vector<double> &t, 
//                   const std::vector<double> &y, 
//                   const int method, 
//                   const std::vector<double> &boundary_condition)
//     : _t(t)
//     , _y(y)
//     , _method(method)
//     , _boundary_condition(boundary_condition)
// {
//     if (!std::is_sorted(t.begin(), t.end()))
//     {
//         std::vector<int> idx(t.size());
//         for (size_t i=0; i<t.size(); ++i)
//             idx[i] = i;
//         std::sort(idx.begin(), idx.end(), [&](int i, int j){return t[i] < t[j];});
//         std::vector<double> t_sorted(t.size()), y_sorted(t.size());
//         for (size_t i=0; i<t.size(); ++i)
//         {
//             t_sorted[i] = t[idx[i]];
//             y_sorted[i] = y[idx[i]];
//         }
//         _t = t_sorted;
//         _y = y_sorted;
//     }
//     interpolate(t, y, method, boundary_condition); // interpolate
// } 

template <>
void 
PPInterpolate<3>::interpolate(
                  const std::vector<double> &t, // nodes
                  const std::vector<double> &y, // values
                  const int method, // 0 for periodic, 1 for complete, 2 for natural, 3 for not-a-knot.
                  const std::vector<double> &boundary_condition) // boundary condition
{
    if (method == 0)
    {
        // To be done
    }
    else if (method == 1)
    {
        int t_size = t.size();
        std::vector<double> K_i(t_size-1, 0);
        for (int i=0; i<t_size-1; ++i)
        {
            K_i[i] = (y[i+1]-y[i])/(t[i+1]-t[i]);
        }

        Eigen::VectorXd b = Eigen::VectorXd::Zero(t_size - 2);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(t_size - 2, t_size - 2);
        for (int i = 0; i < t_size - 2; ++i)
        {
            A(i,i) = 2;
            double mu_i = (t[i+1] - t[i]) / (t[i+2] - t[i]);
            double lambda_i = (t[i+2] - t[i+1]) / (t[i+2] - t[i]);
            b(i) = 3*mu_i*K_i[i+1] + 3*lambda_i*K_i[i];
            if (i == 0)
            {
                A(i, i+1) = mu_i;
            }
            else if (i < t_size - 3)
            {
                A(i, i-1) = lambda_i;
                A(i, i+1) = mu_i;
            }
            else 
            {
                A(i, i-1) = lambda_i;
            }
        }

        // solve Ax = b
        Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);


        std::vector<std::vector<double>> coeffs(t_size-1, std::vector<double>(4,0));
        for (int i=0; i<t_size-1; ++i)
        {
            coeffs[i][0] = y[i];
            if (i == 0)
            {
                double tmp = t[i+1]-t[i];
                coeffs[i][1] = boundary_condition[0];
                coeffs[i][2] = (3*K_i[i]-2*boundary_condition[0]-c(0))/tmp;
                coeffs[i][3] = (boundary_condition[0]-2*K_i[i]+c(0))/(tmp*tmp);
            }
            else if (i < t_size-2)
            {
                double tmp = t[i+1]-t[i];
                coeffs[i][1] = c(i-1);
                coeffs[i][2] = (3*K_i[i]-2*c(i-1)-c(i))/tmp;
                coeffs[i][3] = (c(i-1)-2*K_i[i]+c(i))/(tmp*tmp);
            }
            else
            {
                double tmp = t[i+1]-t[i];
                coeffs[i][1] = c(i-1);
                coeffs[i][2] = (3*K_i[i]-2*c(i-1)-boundary_condition[1])/tmp;
                coeffs[i][3] = (c(i-1)-2*K_i[i]+boundary_condition[1])/(tmp*tmp);
            }
        }
        PPoly tmp = PPoly(coeffs, t);
        poly = tmp;

        // return;
    }
    else if (method == 2)
    {
        /**
         * Natural Cubic Spline
         */
        int t_size = t.size();
        std::vector<double> K_i(t_size-1, 0);
        for (int i=0; i<t_size-1; ++i)
        {
            K_i[i] = (y[i+1]-y[i])/(t[i+1]-t[i]);
        }
        std::vector<double> J_i(t_size-2, 0);
        for (int i=0; i<t_size-2; ++i)
        {
            J_i[i] = (K_i[i+1]-K_i[i])/(t[i+2]-t[i]);
        }

        Eigen::VectorXd b = Eigen::VectorXd::Zero(t_size - 2);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(t_size - 2, t_size - 2);
        for (int i = 0; i < t_size - 2; ++i)
        {
            b(i) = 6*J_i[i];

            double mu_i = (t[i+1] - t[i]) / (t[i+2] - t[i]);
            double lambda_i = (t[i+2] - t[i+1]) / (t[i+2] - t[i]);
            A(i,i) = 2;
            if (i == 0)
            {
                A(i, i+1) = lambda_i;
            }
            else if (i < t_size - 3)
            {
                A(i, i-1) = mu_i;
                A(i, i+1) = lambda_i;
            }
            else
            {
                A(i, i-1) = mu_i;
            }
        }

        // solve Ax = b
        Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

        std::vector<std::vector<double>> coeffs(t.size()-1, std::vector<double>(4,0));
        for (int i=0; i<t_size-1; ++i)
        {
            coeffs[i][0] = y[i];
            if (i == 0)
            {
                double tmp = t[i+1]-t[i];
                coeffs[i][1] = K_i[i]-1.0/6.0 * c(0) * tmp;
                coeffs[i][2] = 0;
                coeffs[i][3] = c(0) / tmp;
            }
            else if (i < t_size-2)
            {
                double tmp = t[i+1]-t[i];
                coeffs[i][1] = K_i[i] - 1.0/6.0 * (c(i+1) + 2*c(i)) * tmp;
                coeffs[i][2] = c(i);
                coeffs[i][3] = (c(i+1) - c(i)) / tmp;
            }
            else
            {
                double tmp = t[i+1]-t[i];
                coeffs[i][1] = K_i[i] - 1.0/3.0 * c(i) * tmp;
                coeffs[i][2] = 0;
                coeffs[i][3] = -c(i) / tmp;
            }
        }

        poly = PPoly(coeffs, t); // interpolate
    }
    // return;
}

// template <>
// double
// PPInterpolate<1>::operator()(const double x) const
// {
//     return poly(x);
// }

// template <>
// double
// PPInterpolate<2>::operator()(const double x) const
// {
//     return poly(x);
// }

// template <>
// double
// PPInterpolate<3>::operator()(const double x) const
// {
//     return poly(x);
// }



template class PPInterpolate<1>;
template class PPInterpolate<2>;
template class PPInterpolate<3>;