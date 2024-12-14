#include "BInterpolate.hpp"
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
    std::cout << "======================C======================" << std::endl;
    std::vector<double> boundary = {10.0 / (26 * 26), -10.0 / (26 * 26)};
    std::vector<double> t1(11);
    std::vector<double> y1(11);
    for (int i = 0; i < 11; i++)
    {
        t1[i] = -5.0 + i;
        y1[i] = f(t1[i]);
    }

    BInterpolate<3> b1(t1, y1, 1, boundary);
    std::string     filename = "./result/C_" + std::to_string(3) + ".txt";
    std::ofstream   ofs(filename);
    for (double x = -5.0; x <= 5.0; x += 0.01)
    {
        ofs << x << "," << b1(x) << std::endl;
    }
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

    filename = "./result/C_" + std::to_string(2) + ".txt";
    ofs      = std::ofstream(filename);
    for (double x = -5.0; x <= 5.0; x += 0.01)
    {
        ofs << x << "," << b2(x) << std::endl;
    }
    ofs.close();
    std::cout << "======================C======================" << std::endl;
    return 0;
}