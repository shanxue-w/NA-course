#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"
#include <fstream>
#include <iostream>
#include <string>

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
    PPInterpolate<3, double> px(t, y, 2);

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

    fout << "max_diff at " << t_mid[max_diff_index] << " is " << diff[max_diff_index] << std::endl;
    fout.close();
    std::cout << "max_diff at " << t_mid[max_diff_index] << " is " << diff[max_diff_index] << std::endl;

    BInterpolate<3, double> bx(t, y, 2);
    filename = "result/A_" + std::to_string(N) + "_bspline.txt";
    fout.open(filename);
    for (double x = -1.0; x <= 1.0; x += 0.01)
    {
        fout << x << "," << bx(x) << std::endl;
    }
    fout.close();
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
    std::cout << "======================A======================" << std::endl;
    return 0;
}