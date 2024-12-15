#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"
#include <fstream>
double
func1(double x)
{
    // 3+8x+4x^2+6x^3
    return 3 + x * (8 + x * (4 + 6 * x));
}

double
func2(double x)
{
    // 1/(1 + e^(-x))
    return 1 / (1 + exp(-x));
}

double
func3(double x)
{
    // sin(x) + e^(-x)
    return sin(x) + exp(-x);
}

int
main()
{
    int                 N = 10;
    std::vector<double> t(N), y(N);
    std::vector<double> boundary_condition = {0, 0};

    std::string   filename = "./result/SplineTest1.txt";
    std::ofstream fout(filename);
    // test function 1
    for (int i = 0; i < N; ++i)
    {
        t[i] = i * 0.1;
        y[i] = func1(t[i]);
    }
    boundary_condition[0] = 8.0;
    boundary_condition[1] = 8.0 + 8.0 * t[N - 1] + 24.0 * t[N - 1] * t[N - 1];

    PPInterpolate<3, double> pp1(t, y, 1, boundary_condition);
    BInterpolate<3, double>  bs1(t, y, 1, boundary_condition);
    for (double i = 0; i < 0.1 * (N - 1); i += 1e-2)
    {
        // std::cout << t[i] << " " << y[i] << " " << pp1(t[i]) << " " << bs1(t[i]) << std::endl;
        fout << i << "," << func1(i) << "," << pp1(i) << "," << bs1(i) << std::endl;
    }
    fout.close();

    filename = "./result/SplineTest2.txt";
    fout.open(filename);
    for (int i = 0; i < N; i++)
    {
        t[i] = -5 + 10.0 / (N - 1) * i;
        y[i] = func2(t[i]);
    }
    boundary_condition[0] = 1.0 / (exp(-5) + exp(5) + 2);
    boundary_condition[1] = 1.0 / (exp(-5) + exp(5) + 2);
    PPInterpolate<3, double> pp2(t, y, 1, boundary_condition);
    BInterpolate<3, double>  bs2(t, y, 1, boundary_condition);
    for (double i = -5; i < 5; i += 1e-2)
    {
        fout << i << "," << func2(i) << "," << pp2(i) << "," << bs2(i) << std::endl;
    }
    fout.close();

    filename = "./result/SplineTest3.txt";
    fout.open(filename);
    for (int i = 0; i < N; i++)
    {
        t[i] = -2 * M_PI + 4 * M_PI / (N - 1) * i;
        y[i] = func3(t[i]);
    }
    boundary_condition[0] = cos(-2 * M_PI) - exp(2 * M_PI);
    boundary_condition[1] = cos(2 * M_PI) - exp(-2 * M_PI);
    PPInterpolate<3, double> pp3(t, y, 1, boundary_condition);
    BInterpolate<3, double>  bs3(t, y, 1, boundary_condition);
    for (double i = -2 * M_PI; i < 2 * M_PI; i += 1e-2)
    {
        fout << i << "," << func3(i) << "," << pp3(i) << "," << bs3(i) << std::endl;
    }
    fout.close();

    return 0;
}