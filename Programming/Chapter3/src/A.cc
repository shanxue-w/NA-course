#include "BInterpolate.hpp"
#include "MinLeastSquare.hpp"
#include "PPInterpolate.hpp"
#include <fstream>
#include <iostream>
#include <string>

std::vector<double> x_lists   = {std::log(6), std::log(11), std::log(21), std::log(41), std::log(81)};
std::vector<double> y_lists   = {0.0, 0.0, 0.0, 0.0, 0.0};
int                 index_ofy = 0;

double
f(double x)
{
    return 1.0 / (1.0 + 25.0 * x * x);
}

void
run(int N)
{
    std::vector<double> t(N);
    std::vector<double> y(N);
    for (int i = 0; i < N; i++)
    {
        t[i] = -1.0 + 2.0 * i / (N - 1);
        y[i] = f(t[i]);
    }
    // use natural cubic spline interpolation
    PPInterpolate<3> px(t, y, 2);

    std::string   filename = "result/A_" + std::to_string(N) + ".txt";
    std::ofstream fout(filename);

    for (double x = -1.0; x <= 1.0; x += 0.01)
    {
        fout << x << "," << px(x) << std::endl;
    }
    fout.close();
    filename = "result/A_" + std::to_string(N) + "_mid.txt";
    fout.open(filename);
    std::vector<double> diff(N - 1);
    std::vector<double> t_mid(N - 1);
    for (int i = 0; i < N - 1; i++)
    {
        double mid = (t[i] + t[i + 1]) / 2;
        // std::cout << mid << "," << px(mid) << "," << f(mid) << std::endl;
        fout << mid << "," << px(mid) << "," << f(mid) << std::endl;
        // std::cout << "diff: " << px(mid) - f(mid) << std::endl;
        fout << "diff: " << px(mid) - f(mid) << std::endl;
        diff[i]  = std::abs(px(mid) - f(mid));
        t_mid[i] = mid;
    }
    // find the index of the maximum difference
    int max_diff_index = std::max_element(diff.begin(), diff.end()) - diff.begin();
    y_lists[index_ofy] = std::log(diff[max_diff_index]);
    index_ofy++;

    fout << "max_diff at " << t_mid[max_diff_index] << " is " << diff[max_diff_index] << std::endl;
    fout.close();
    std::cout << "max_diff at " << t_mid[max_diff_index] << " is " << diff[max_diff_index] << std::endl;
}

int
main(void)
{
    std::cout << "======================A======================" << std::endl;
    run(6);
    run(11);
    run(21);
    run(41);
    run(81);
    double k = MinLeastSquare(x_lists, y_lists);
    std::cout << "k = " << k << std::endl;
    std::cout << "======================A======================" << std::endl;
    return 0;
}