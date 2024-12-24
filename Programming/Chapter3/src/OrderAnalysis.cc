#include "BInterpolate.hpp"
#include "MinLeastSquare.hpp"
#include "PPInterpolate.hpp"

double
f(double x)
{
    return 1.0 / (1.0 + 25.0 * x * x);
}

int
main(void)
{
    std::cout << "======================OrderAnalysis======================" << std::endl;
    int n = 30;

    std::vector<int>    N_lists;
    std::vector<double> x_lists;
    for (int i = 0; i < n; i++)
    {
        N_lists.push_back(10 * i + 10);
        x_lists.push_back(std::log(10 * i + 9));
    }
    std::vector<double> boundary_condition = {50.0 / (26.0 * 26.0), -50.0 / (26.0 * 26.0)};
    std::vector<double> y_lists_order1(n, 0.0);
    std::vector<double> y_lists_order2(n, 0.0);
    std::vector<double> y_lists_order3(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        std::vector<double> x(N_lists[i], 0.0);
        std::vector<double> y(N_lists[i], 0.0);
        for (int j = 0; j < N_lists[i]; j++)
        {
            x[j] = -1.0 + 2.0 * j / (N_lists[i] - 1);
            y[j] = f(x[j]);
        }

        PPInterpolate<1> px1(x, y, 1, boundary_condition);
        PPInterpolate<2> px2(x, y, 1, boundary_condition);
        PPInterpolate<3> px3(x, y, 1, boundary_condition);

        double max_error1 = 0.0, max_error2 = 0.0, max_error3 = 0.0;
        for (double tmp = -1.0; tmp <= 1.0; tmp += 0.01)
        {
            max_error1 = std::max(max_error1, std::abs(px1(tmp) - f(tmp)));
            max_error2 = std::max(max_error2, std::abs(px2(tmp) - f(tmp)));
            max_error3 = std::max(max_error3, std::abs(px3(tmp) - f(tmp)));
        }

        y_lists_order1[i] = std::log(max_error1);
        y_lists_order2[i] = std::log(max_error2);
        y_lists_order3[i] = std::log(max_error3);
    }
    // y = A h^j = A (2/(N-1)) ^j
    // log(y) = C - j log(N-1)
    std::cout << "Order 1 convergence rate: " << -MinLeastSquare(x_lists, y_lists_order1) << std::endl;
    std::cout << "Order 2 convergence rate: " << -MinLeastSquare(x_lists, y_lists_order2) << std::endl;
    std::cout << "Order 3 convergence rate: " << -MinLeastSquare(x_lists, y_lists_order3) << std::endl;

    std::cout << "======================OrderAnalysis======================" << std::endl;
    return 0;
}