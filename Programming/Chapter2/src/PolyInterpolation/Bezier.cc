/**
 * @file Bezier.cc
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief Implementation of Bezier Curve.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "Bezier.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

/**
 * @brief Construct a new Bezier:: Bezier object
 * 
 * To use the Bezier Curve, we need the following data:
 * \f{eqnarray*}{
 * \mathbf{P}(x_j) = (x_j, y_j), j = 0, 1, 2, \cdots, m-1 \\
 * \mathbf{P}^{\prime}(x_j) = (P^{\prime}_x(x_j), P^{\prime}_y(x_j)), j = 0, 1, 2, \cdots, m-1
 * \f}
 * 
 * @param xData the x values of the data points
 * @param yData the y values of the data points
 * @param PxData the x values of the derivative of the data points, can multiply by the same scalar with PyData.
 * @param PyData the y values of the derivative of the data points, can multiply by the same scalar with PxData.
 * 
 * @code {.cc}
 * std::vector<double> xData = {0, 1, 2, 3};
 * std::vector<double> yData = {0, 1, 4, 9};
 * std::vector<double> PxData = {1, 2, 3, 4};
 * std::vector<double> PyData = {1, 2, 3, 4};
 * Bezier Bezier(xData, yData, PxData, PyData);
 * @endcode
 */
Bezier::Bezier(
    const std::vector<double> &xData, 
    const std::vector<double> &yData,
    const std::vector<double> &PxData,
    const std::vector<double> &PyData
    )
{
    m = xData.size();
    for (int i = 0; i < m; i++)
    {
        x_b0_lists.push_back(xData[i]);
        x_b1_lists.push_back(xData[i] + 1.0/3.0 * PxData[i]);
        x_b2_lists.push_back(xData[(i+1)%m] - 1.0/3.0 * PxData[(i+1)%m]);
        x_b3_lists.push_back(xData[(i+1)%m]);

        y_b0_lists.push_back(yData[i]);
        y_b1_lists.push_back(yData[i] + 1.0/3.0 * PyData[i]);
        y_b2_lists.push_back(yData[(i+1)%m] - 1.0/3.0 * PyData[(i+1)%m]);
        y_b3_lists.push_back(yData[(i+1)%m]);
    }
}

/**
 * @brief Calculate the value of the Bezier Curve at the point i. 
 * 
 * \f[
 * x = x(t), y = y(t), t \in [i, i+1]
 * \f]
 * 
 * @param i the value of t
 * @return std::vector<double> the value of the Bezier Curve at the point i
 * 
 * @code {.cc}
 * Bezier Bezier;
 * std::vector<double> point = Bezier(1.5);
 * @endcode
 */
std::vector<double> Bezier::operator()(double i) const
{
    int k = std::floor(i);
    double tmp = i - k;
    k = k % m;
    double x = x_b0_lists[k] * B0(tmp) + x_b1_lists[k] * B1(tmp) + x_b2_lists[k] * B2(tmp) + x_b3_lists[k] * B3(tmp);
    double y = y_b0_lists[k] * B0(tmp) + y_b1_lists[k] * B1(tmp) + y_b2_lists[k] * B2(tmp) + y_b3_lists[k] * B3(tmp);
    std::vector<double> point = {x, y};
    return point; 
}

