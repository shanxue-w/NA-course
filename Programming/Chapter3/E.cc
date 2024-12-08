#include "BSpline.hpp"
#include "Curve.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/**
 * @details
 *
 * For \f$r(t) = (\sin t+t\cos t, \cos t - t\sin t)\f$,
 * we have \f$ r^{\prime}(t) = (2\cos t - t\sin t, -2\sin t - t\cos t) \f$.
 * So, \f$ r^{\prime}(t) \cdot r^{\prime}(t) = 4 + t^2 \f$.
 * So, \f$ s(t) = \int \sqrt{4+t^2}\mathrm{d}t \f$.
 * We can easily get
 * \f$ s(t) = \frac{1}{2}t\sqrt{4+t^2} + 2\sinh^{-1}\frac{t}{2} \f$.
 * \f$ \sinh^{-1}(x) = \log(x + \sqrt{x^2+1}) \f$.
 *
 */

void
runA(int N)
{
    std::vector<double> t(N + 1);
    // std::vector<double> x(N + 1), y(N + 1);
    std::vector<std::vector<double>> xy(2, std::vector<double>(N + 1));
    const double                     sqrt3 = std::sqrt(3.0);
    for (int i = 0; i < N + 1; ++i)
    {
        t[i]     = i * 2 * M_PI / N;
        xy[0][i] = sqrt3 * sin(t[i]);
        xy[1][i] =
            2.0 / 3.0
            * (sqrt3 * cos(t[i]) + std::sqrt(sqrt3 * std::abs(sin(t[i]))));
    }

    Curve<3>      curve(t, xy);
    std::string   filename = "./result/E_heart_B" + std::to_string(N) + ".txt";
    std::ofstream fout(filename);
    auto          bspline = curve.BSpline2d();
    BSpline       bx      = bspline[0];
    BSpline       by      = bspline[1];
    for (double x = 0.0; x <= 2 * M_PI; x += 1e-4)
    {
        fout << bx(x) << "," << by(x) << std::endl;
    }
    fout.close();

    auto  ppoly = curve.PPoly2d();
    PPoly px    = ppoly[0];
    PPoly py    = ppoly[1];
    filename    = "./result/E_heart_PP" + std::to_string(N) + ".txt";
    fout.open(filename);

    for (double x = 0.0; x <= 2 * M_PI; x += 1e-4)
    {
        fout << px(x) << "," << py(x) << std::endl;
    }
    fout.close();
    return;
}

void
runB(int N)
{
    std::vector<double>              t(N + 1);
    std::vector<std::vector<double>> xy(2, std::vector<double>(N + 1));
    for (int i = 0; i < N + 1; ++i)
    {
        t[i]     = i * 2 * M_PI / N;
        xy[0][i] = sin(t[i]) + t[i] * cos(t[i]);
        xy[1][i] = cos(t[i]) - t[i] * sin(t[i]);
    }

    Curve<3>      curve(t, xy);
    std::string   filename = "./result/E_2_B" + std::to_string(N) + ".txt";
    std::ofstream fout(filename);
    auto          bspline = curve.BSpline2d();
    BSpline       bx      = bspline[0];
    BSpline       by      = bspline[1];
    for (double x = 0.0; x <= 2 * M_PI; x += 0.001)
    {
        fout << bx(x) << "," << by(x) << std::endl;
    }
    fout.close();

    auto  ppoly = curve.PPoly2d();
    PPoly px    = ppoly[0];
    PPoly py    = ppoly[1];
    filename    = "./result/E_2_PP" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = 0.0; x <= 2 * M_PI; x += 0.001)
    {
        fout << px(x) << "," << py(x) << std::endl;
    }
    fout.close();
    return;
}

void
runC(int N)
{
    std::vector<double>              t(N + 1);
    std::vector<std::vector<double>> xy(3, std::vector<double>(N + 1));

    for (int i = 0; i < N + 1; ++i)
    {
        t[i]     = i * 2 * M_PI / N;
        xy[0][i] = sin(cos(t[i])) * cos(sin(t[i]));
        xy[1][i] = sin(cos(t[i])) * sin(sin(t[i]));
        xy[2][i] = cos(cos(t[i]));
    }

    Curve<3>      curve(t, xy);
    std::string   filename = "./result/E_3_B" + std::to_string(N) + ".txt";
    std::ofstream fout(filename);
    auto          bspline = curve.BSpline3d();
    BSpline       bx      = bspline[0];
    BSpline       by      = bspline[1];
    BSpline       bz      = bspline[2];
    for (double x = 0.0; x <= 2 * M_PI; x += 0.001)
    {
        fout << bx(x) << "," << by(x) << "," << bz(x) << std::endl;
    }
    fout.close();

    auto  ppoly = curve.PPoly3d();
    PPoly px    = ppoly[0];
    PPoly py    = ppoly[1];
    PPoly pz    = ppoly[2];
    filename    = "./result/E_3_PP" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = 0.0; x <= 2 * M_PI; x += 0.001)
    {
        fout << px(x) << "," << py(x) << "," << pz(x) << std::endl;
    }
    fout.close();

    auto bsplineball = curve.BSplineBall();
    filename         = "./result/E_3_Ball" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = 0.0; x <= 2 * M_PI; x += 0.001)
    {
        std::vector<double> p = bsplineball(x);
        fout << p[0] << "," << p[1] << "," << p[2] << std::endl;
    }
    fout.close();

    auto ppolyball = curve.PPolyBall();
    filename       = "./result/E_3_PPBall" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = 0.0; x <= 2 * M_PI; x += 0.001)
    {
        std::vector<double> p = ppolyball(x);
        fout << p[0] << "," << p[1] << "," << p[2] << std::endl;
    }
    fout.close();
    return;
}

int
main(void)
{
    /**
     * @details
     * For the heart curve, we know the equation is
     * \f$ x^2 + (\frac{3}{2}y - \sqrt{|x|})^2 = 3. \f$
     * then we can get
     * \f[
     * x = \sqrt{3} \sin t, \quad
     * y = \frac{2}{3} (\sqrt{3} \cos t + \sqrt{\sqrt{3} |\sin t|})
     * \f].
     */
    runA(10);
    runA(40);
    runA(160);
    runA(500);
    runB(100);
    runC(100);
}