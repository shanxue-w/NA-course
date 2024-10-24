#include <iostream>
#include "Bezier.hpp"
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

int m = 10;
double sqrt3 = std::sqrt(3.0);
std::vector<double> xData(m);
std::vector<double> yData(m);
std::vector<double> PxData(m);
std::vector<double> PyData(m);

double get_y_positive(double x)
{
    return 2.0/3.0 * (std::sqrt(std::abs(x)) + std::sqrt(3.0-x*x));
}

double get_y_negative(double x)
{
    return 2.0/3.0 * (std::sqrt(std::abs(x)) - std::sqrt(3.0-x*x));
}

double prime(double x, double y)
{
    if (x > 0)
    {
        double sqrtx = std::sqrt(x);
        return 1.0/3.0 * (1/sqrtx - 2*x/(1.5*y-sqrtx));
    }
    else if (x < 0)
    {
        double sqrtx = std::sqrt(-x);
        return -1.0/3.0 * (1/sqrtx + 2*x/(1.5*y-sqrtx));
    }
    return 0;
}

void setup(int m)
{
    xData.resize(m);
    yData.resize(m);
    PxData.resize(m);
    PyData.resize(m);
    double delta_x = (4.0 * sqrt3- 4e-3) / (m-2);
    for (int i=0; i< m/2; i++)
    {
        double x = -sqrt3 + 1e-3 + i * delta_x;
        xData[i] = x;
        yData[i] = get_y_positive(x);
        PxData[i] = delta_x / 15.0;
        PyData[i] = prime(x, yData[i]) * delta_x / 15.0;
    }
    for (int i=0; i< m/2; i++)
    {
        double x = sqrt3 - 1e-3 - i * delta_x;
        xData[i+m/2] = x;
        yData[i+m/2] = get_y_negative(x);
        PxData[i+m/2] = -delta_x / 15.0;
        PyData[i+m/2] = -prime(x, yData[i+m/2]) * delta_x / 15.0;
    }
}

int main(void)
{
    setup(10);
    Bezier Bezier1(xData, yData, PxData, PyData);
    std::string filename = "./data/F_Bezier_" + std::to_string(10) + ".txt";
    std::ofstream output(filename);
    for (double i = 0; i <= 10; i+=0.01)
    {
        std::vector<double> point = Bezier1(i);
        output << point[0] << "," << point[1] << std::endl;
    }
    output.close();

    setup(40);
    Bezier Bezier2(xData, yData, PxData, PyData);
    filename = "./data/F_Bezier_" + std::to_string(40) + ".txt";
    std::ofstream output2(filename);
    for (double i = 0; i <= 40; i+=0.01)
    {
        std::vector<double> point = Bezier2(i);
        output2 << point[0] << "," << point[1] << std::endl;
    }
    output2.close();

    setup(160);
    Bezier Bezier3(xData, yData, PxData, PyData);
    filename = "./data/F_Bezier_" + std::to_string(160) + ".txt";
    std::ofstream output3(filename);
    for (double i = 0; i <= 160; i+=0.01)
    {
        std::vector<double> point = Bezier3(i);
        output3 << point[0] << "," << point[1] << std::endl;
    }
    output3.close();
    return 0;
}