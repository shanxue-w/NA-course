#include "BInterpolate.hpp"
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

double
f(double x)
{
    return 1.0 / (1.0 + x * x);
}

int
main(void)
{
    std::cout << "======================D======================" << std::endl;
    std::vector<double> boundary = {10.0 / (26 * 26), -10.0 / (26 * 26)};
    std::vector<double> t1(11);
    std::vector<double> y1(11);
    for (int i = 0; i < 11; i++)
    {
        t1[i] = -5.0 + i;
        y1[i] = f(t1[i]);
    }

    BInterpolate<3> b1(t1, y1, 1, boundary);
    std::string             filename = "./result/D_" + std::to_string(3) + ".txt";
    std::ofstream           ofs(filename);
    std::vector<double>     x_lists = {-3.5, -2.5, -1.5, -0.5, 0.0, 0.5, 1.5, 2.5, 3.5};
    for (double x : x_lists)
    {
        ofs << x << ", " << std::abs(b1(x) - f(x)) << std::endl;
        std::cout << x << ", " << std::abs(b1(x) - f(x)) << std::endl;
    }
    std::cout << std::endl;
    ofs.close();

    std::vector<double> t2(10);
    std::vector<double> y2(10);
    for (int i = 0; i < 10; i++)
    {
        t2[i] = -4.5 + i;
        y2[i] = f(t2[i]);
    }
    boundary = {f(-5.0), f(5.0)};

    BInterpolate<2> b2(t2, y2, 2, boundary);

    filename = "./result/D_" + std::to_string(2) + ".txt";
    ofs      = std::ofstream(filename);

    for (double x : x_lists)
    {
        ofs << x << ", " << std::abs(b2(x) - f(x)) << std::endl;
        std::cout << x << ", " << std::abs(b2(x) - f(x)) << std::endl;
    }
    ofs.close();

    double abs_max1 = 0.0;
    double abs_max2 = 0.0;
    for (double x = -5.0; x <= 5.0; x += 0.01)
    {
        abs_max1 = std::max(abs_max1, std::abs(b1(x) - f(x)));
        abs_max2 = std::max(abs_max2, std::abs(b2(x) - f(x)));
    }
    filename = "./result/D_abs_max.txt";
    ofs      = std::ofstream(filename);
    ofs << "max error for degree 3: " << abs_max1 << std::endl;
    std::cout << "max error for degree 3: " << abs_max1 << std::endl;
    ofs << "max error for degree 2: " << abs_max2 << std::endl;
    std::cout << "max error for degree 2: " << abs_max2 << std::endl;
    ofs.close();
    std::cout << "======================D======================" << std::endl;
    return 0;
}