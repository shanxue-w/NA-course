// 高阶PP-form和 高阶B-form的测试代码
// n阶样条至少应该能够完美拟合n-1阶多项式，这是一个测试样例
// 其它一些简单曲线的拟合样例
// 高阶只实现了两个，0-周期样条，1-完全样条（所有边界条件在起始点）
// 超出区间范围是强制把样条值设为0，测试时注意。

#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"
#include <cmath>

int
main(void)
{
    std::cout << "======================HIGH======================" << std::endl;
    int    N        = 50;
    double maxerror = 0.0;

    // f(x) = 1 + x + x^2 + x^3
    auto                function1 = [](double x) { return 1 + x * (1 + x * (1 + x)); };
    std::vector<double> x(N), y(N);
    std::vector<double> boundary_condition(10, 0.0);
    for (int i = 0; i < N; i++)
    {
        x[i] = i;
        y[i] = function1(i);
    }
    boundary_condition[0] = 1.0;
    boundary_condition[1] = 2.0;
    boundary_condition[2] = 6.0;
    BInterpolate<4, double>  b1(x, y, 1, boundary_condition);
    PPInterpolate<4, double> p1(x, y, 1, boundary_condition);
    maxerror = 0.0;
    for (double x = 0.0; x < N - 1; x += 1e-2)
    {
        double error = std::abs(b1(x) - function1(x));
        if (error > maxerror)
            maxerror = error;
    }
    std::cout << "max error of BInterpolate<4, double>: " << maxerror << std::endl;
    maxerror = 0.0;
    for (double x = 0.0; x < N - 1; x += 1e-2)
    {
        double error = std::abs(p1(x) - function1(x));
        if (error > maxerror)
            maxerror = error;
    }
    std::cout << "max error of PPInterpolate<4, double>: " << maxerror << std::endl;

    // f(x) = sin(x)
    mpf_set_default_prec(1000);
    for (int i = 0; i < N; i++)
    {
        x[i] = 2.0 * M_PI * i / (N - 1);
        y[i] = sin(x[i]);
    }
    BInterpolate<5, double> b2(x, y, 0);
    maxerror = 0.0;
    for (double x = 0.0; x < 2.0 * M_PI; x += 1e-2)
    {
        double error = std::abs(b2(x) - sin(x));
        if (error > maxerror)
            maxerror = error;
    }
    std::cout << "sin(x) max error of BInterpolate<5, double>: " << maxerror << std::endl;
    std::vector<mpf_class> x2(N), y2(N);
    for (int i = 0; i < N; i++)
    {
        x2[i] = 2.0 * M_PI * i / (N - 1);
        y2[i] = sin(x2[i].get_d());
    }
    PPInterpolate<5, mpf_class> p2(x2, y2, 0);
    mpf_class                   max_err = mpf_class(0.0);
    max_err                             = mpf_class(0.0);
    for (double x = 0.0; x < 2.0 * M_PI; x += 1e-2)
    {
        mpf_class error = abs(p2(x) - sin(x));
        if (error > max_err)
            max_err = error;
    }
    std::cout << "sin(x) max error of PPInterpolate<5, mpf_class>: " << max_err << std::endl;

    // f(x) = sin(x) + cos(x)
    std::vector<double> x4(N), y4(N);
    for (int i = 0; i < N; i++)
    {
        x4[i] = i * 2.0 * M_PI / (N - 1);
        y4[i] = sin(x4[i]) + cos(x4[i]);
    }
    BInterpolate<7, double> b4(x4, y4, 0);
    maxerror = 0.0;
    for (double x = 0.0; x < 2.0 * M_PI; x += 1e-2)
    {
        double error = std::abs(b4(x) - (sin(x) + cos(x)));
        if (error > maxerror)
            maxerror = error;
    }
    std::cout << "sin(x)+cos(x) max error of BInterpolate<7, double>: " << maxerror << std::endl;

    mpf_set_default_prec(1000);
    std::vector<mpf_class> x3(N), y3(N);
    for (int i = 0; i < N; i++)
    {
        x3[i] = i * 2.0 * M_PI / (N - 1);
        y3[i] = sin(x3[i].get_d()) + cos(x3[i].get_d());
    }
    PPInterpolate<7, mpf_class> p3(x3, y3, 0);
    max_err = mpf_class(0.0);
    for (double x = 0.0; x < 2.0 * M_PI; x += 1e-2)
    {
        mpf_class error = abs(p3(x) - (sin(x) + cos(x)));
        if (error > max_err)
            max_err = error;
    }
    std::cout << "sin(x)+cos(x) max error of PPInterpolate<7, mpf_class>: " << max_err << std::endl;
    std::cout << "======================HIGH======================" << std::endl;
    return 0;
}