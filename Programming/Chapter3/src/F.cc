#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::vector<std::vector<double>>
dividedDifference(std::vector<double> &t, double &x)
{
    int                              t_size = t.size();
    int                              n      = t_size - 2;
    std::vector<std::vector<double>> dd(t_size, std::vector<double>({}));
    double                           tmp = t[t_size - 1] - t[0];
    for (int i = 0; i < t_size; i++)
    {
        // compute (t[i]-x)_{+}^n
        if (t[i] >= x)
        {
            dd[0].push_back(std::pow((t[i] - x), n) * tmp);
        }
        else
        {
            dd[0].push_back(0.0);
        }
    }
    // compute divided difference
    for (int j = 1; j < t_size; j++)
    {
        for (int i = 0; i < t_size - j; i++)
        {
            dd[j].push_back((dd[j - 1][i + 1] - dd[j - 1][i]) / (t[i + j] - t[i]));
        }
    }
    return dd;
}

int
main()
{
    std::cout << "======================F======================" << std::endl;
    std::string   filename = "./result/F_1.txt";
    std::ofstream file(filename);

    std::vector<double> t = {0.0, 1.0, 2.0};
    for (double x = -0.5; x < 2.5; x += 0.01)
    {
        auto data = dividedDifference(t, x);
        file << x << "," << data[0][0] << "," << data[0][1] << "," << data[0][2] << "," << data[1][0] << ","
             << data[1][1] << "," << data[2][0] << std::endl;
    }
    file.close();

    filename = "./result/F_2.txt";
    file.open(filename);
    t = {0.0, 1.0, 2.0, 3.0};
    for (double x = -0.5; x < 3.5; x += 0.01)
    {
        auto data = dividedDifference(t, x);
        file << x << "," << data[0][0] << "," << data[0][1] << "," << data[0][2] << "," << data[0][3] << ","
             << data[1][0] << "," << data[1][1] << "," << data[1][2] << "," << data[2][0] << "," << data[2][1] << ","
             << data[3][0] << std::endl;
    }
    file.close();
    std::cout << "======================F======================" << std::endl;
    return 0;
}