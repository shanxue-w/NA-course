/**
 * @file Bezier.hpp
 * @author WangHao (3220104819@zju.edu.cn)
 * @brief The Bezier class, providing the basic operations of Bezier.
 * @version 0.1
 * @date 2024-10-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __Bezier_HPP__
#define __Bezier_HPP__

#include <vector>

class Bezier
{
public:
    Bezier( const std::vector<double> &xData, 
            const std::vector<double> &yData,
            const std::vector<double> &PxData,
            const std::vector<double> &PyData
          );
    std::vector<double> operator()(double i) const; // 计算函数值
    double B0(double x) const { return (1-x)*(1-x)*(1-x); }
    double B1(double x) const { return 3*x*(1-x)*(1-x); }
    double B2(double x) const { return 3*x*x*(1-x); }
    double B3(double x) const { return x*x*x; }
private:
    // std::vector<double> x_lists;
    std::vector<double> x_b0_lists;
    std::vector<double> x_b1_lists;
    std::vector<double> x_b2_lists;
    std::vector<double> x_b3_lists;

    std::vector<double> y_b0_lists;
    std::vector<double> y_b1_lists;
    std::vector<double> y_b2_lists;
    std::vector<double> y_b3_lists;
    int m;
};



#endif