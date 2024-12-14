#ifndef __MINLEASTSQUARE_HPP__
#define __MINLEASTSQUARE_HPP__

#include <Eigen/Dense>
#include <cmath>
#include <vector>
// 最小二乘法求斜率
double
MinLeastSquare(const std::vector<double> &x, const std::vector<double> &y)
{
    Eigen::MatrixXd A(x.size(), 2);
    A.col(0) = Eigen::VectorXd::Ones(x.size());
    A.col(1) = Eigen::VectorXd::Map(&x[0], x.size());
    Eigen::VectorXd b(y.size());
    b                  = Eigen::VectorXd::Map(&y[0], y.size());
    Eigen::VectorXd x_ = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    return x_[1];
}

#endif