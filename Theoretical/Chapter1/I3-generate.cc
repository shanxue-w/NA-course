#include <iostream>
#include <functional>
#include <iomanip>
double f(double x)
{
    return 4*x*x*x - 2*x*x + 3;
}

double df(double x)
{
    return 12*x*x - 4*x;
}

int main(void)
{
    double x_0 = -1.0;
    // 保留十位小数
    std::cout << 0 << " & " << std::fixed << std::setprecision(9) << x_0 << " & "<<std::fixed << std::setprecision(9) << f(x_0) << " \\\\ \\hline" << std::endl;
    for(auto i = 1; i < 20; i++)
    {
        x_0 = x_0 - f(x_0)/df(x_0);
        std::cout << i << " & " << std::fixed << std::setprecision(9) << x_0 << " & " <<std::setprecision(9) << f(x_0) << " \\\\ \\hline" << std::endl;
    }
    return 0;
}