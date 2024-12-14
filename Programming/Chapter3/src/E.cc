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
    std::vector<double> t_cum(N + 1, 0.0);
    // std::vector<double> x(N + 1), y(N + 1);
    std::vector<std::vector<double>> xy(2, std::vector<double>(N + 1));
    const double                     sqrt3 = std::sqrt(3.0);
    for (int i = 0; i < N + 1; ++i)
    {
        t[i]     = i * 2 * M_PI / N;
        xy[0][i] = sqrt3 * std::sin(t[i]);
        xy[1][i] = 2.0 / 3.0 * (sqrt3 * cos(t[i]) + std::sqrt(sqrt3 * std::abs(std::sin(t[i]))));
    }
    for (int i = 1; i < N + 1; ++i)
    {
        double dx = xy[0][i] - xy[0][i - 1];
        double dy = xy[1][i] - xy[1][i - 1];
        t_cum[i]  = t_cum[i - 1] + std::sqrt(dx * dx + dy * dy);
    }

    Curve<3>      curve(t, xy);
    std::string   filename = "./result/E_heart_B" + std::to_string(N) + ".txt";
    std::ofstream fout(filename);
    auto          bspline = curve.BSpline2d();
    BSpline       bx      = bspline[0];
    BSpline       by      = bspline[1];
    for (double x = 0.0; x <= 2 * M_PI; x += 1e-2)
    {
        fout << bx(x) << "," << by(x) << std::endl;
    }
    fout.close();

    auto  ppoly = curve.PPoly2d();
    PPoly px    = ppoly[0];
    PPoly py    = ppoly[1];
    filename    = "./result/E_heart_PP" + std::to_string(N) + ".txt";
    fout.open(filename);

    for (double x = 0.0; x <= 2 * M_PI; x += 1e-2)
    {
        fout << px(x) << "," << py(x) << std::endl;
    }
    fout.close();

    Curve<3>      curve_cum(t_cum, xy);
    std::string   filename_cum = "./result/E_heart_cum_B" + std::to_string(N) + ".txt";
    std::ofstream fout_cum(filename_cum);
    auto          bspline_cum = curve_cum.BSpline2d();
    BSpline       bx_cum      = bspline_cum[0];
    BSpline       by_cum      = bspline_cum[1];
    for (double x = 0.0; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << bx_cum(x) << "," << by_cum(x) << std::endl;
    }
    fout_cum.close();

    auto  ppoly_cum = curve_cum.PPoly2d();
    PPoly px_cum    = ppoly_cum[0];
    PPoly py_cum    = ppoly_cum[1];
    filename_cum    = "./result/E_heart_cum_PP" + std::to_string(N) + ".txt";
    fout_cum.open(filename_cum);
    for (double x = 0.0; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << px_cum(x) << "," << py_cum(x) << std::endl;
    }
    fout_cum.close();
    return;
}

void
runB(int N)
{
    std::vector<double>              t(N + 1);
    std::vector<double>              t_cum(N + 1, 0.0);
    std::vector<std::vector<double>> xy(2, std::vector<double>(N + 1));
    std::vector<std::vector<double>> boundary_condition(2, std::vector<double>(2));
    boundary_condition[0][0] = 2.0;
    boundary_condition[0][1] = 2.0;
    boundary_condition[1][0] = 0.0;
    boundary_condition[1][1] = -6.0 * M_PI;
    for (int i = 0; i < N + 1; ++i)
    {
        t[i]     = i * 6 * M_PI / N;
        xy[0][i] = sin(t[i]) + t[i] * cos(t[i]);
        xy[1][i] = cos(t[i]) - t[i] * sin(t[i]);
    }
    for (int i = 1; i < N + 1; ++i)
    {
        double dx = xy[0][i] - xy[0][i - 1];
        double dy = xy[1][i] - xy[1][i - 1];
        t_cum[i]  = t_cum[i - 1] + std::sqrt(dx * dx + dy * dy);
    }
    Curve<3>      curve(t, xy, 1, boundary_condition);
    std::string   filename = "./result/E_2_B" + std::to_string(N) + ".txt";
    std::ofstream fout(filename);
    auto          bspline = curve.BSpline2d();
    BSpline       bx      = bspline[0];
    BSpline       by      = bspline[1];
    for (double x = 0.0; x <= 6 * M_PI; x += 1e-2)
    {
        fout << bx(x) << "," << by(x) << std::endl;
    }
    fout.close();

    auto  ppoly = curve.PPoly2d();
    PPoly px    = ppoly[0];
    PPoly py    = ppoly[1];
    filename    = "./result/E_2_PP" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = 0.0; x <= 6 * M_PI; x += 1e-2)
    {
        fout << px(x) << "," << py(x) << std::endl;
    }
    fout.close();

    // ds / dt = sqrt(4 + t^2)
    // so dx/ds = dx/dt * dt/ds = dx/dt * sqrt(4 + t^2)
    boundary_condition[0][0] = 2.0 / 2.0;
    boundary_condition[0][1] = 2.0 / std::sqrt(4 + (6 * M_PI) * (6 * M_PI));
    boundary_condition[1][0] = 0.0 / 2.0;
    boundary_condition[1][1] = -6.0 * M_PI / std::sqrt(4 + (6 * M_PI) * (6 * M_PI));
    Curve<3>      curve_cum(t_cum, xy, 1, boundary_condition);
    std::string   filename_cum = "./result/E_2_cum_B" + std::to_string(N) + ".txt";
    std::ofstream fout_cum(filename_cum);
    auto          bspline_cum = curve_cum.BSpline2d();
    BSpline       bx_cum      = bspline_cum[0];
    BSpline       by_cum      = bspline_cum[1];
    for (double x = 0.0; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << bx_cum(x) << "," << by_cum(x) << std::endl;
    }
    fout_cum.close();

    auto  ppoly_cum = curve_cum.PPoly2d();
    PPoly px_cum    = ppoly_cum[0];
    PPoly py_cum    = ppoly_cum[1];
    filename_cum    = "./result/E_2_cum_PP" + std::to_string(N) + ".txt";
    fout_cum.open(filename_cum);
    for (double x = 0.0; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << px_cum(x) << "," << py_cum(x) << std::endl;
    }
    fout_cum.close();
    return;
}

void
runC(int N)
{
    std::vector<double>              t(N + 1);
    std::vector<double>              t_cum(N + 1, 0.0);
    std::vector<std::vector<double>> xy(3, std::vector<double>(N + 1));

    for (int i = 0; i < N + 1; ++i)
    {
        t[i]     = i * 2 * M_PI / N;
        xy[0][i] = sin(cos(t[i])) * cos(sin(t[i]));
        xy[1][i] = sin(cos(t[i])) * sin(sin(t[i]));
        xy[2][i] = cos(cos(t[i]));
    }
    for (int i = 1; i < N + 1; ++i)
    {
        double dx = xy[0][i] - xy[0][i - 1];
        double dy = xy[1][i] - xy[1][i - 1];
        double dz = xy[2][i] - xy[2][i - 1];
        t_cum[i]  = t_cum[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    Curve<3>      curve(t, xy);
    std::string   filename = "./result/E_3_B" + std::to_string(N) + ".txt";
    std::ofstream fout(filename);
    auto          bspline = curve.BSpline3d();
    BSpline       bx      = bspline[0];
    BSpline       by      = bspline[1];
    BSpline       bz      = bspline[2];
    for (double x = t[0]; x <= t[N]; x += 1e-2)
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
    for (double x = t[0]; x <= t[N]; x += 1e-2)
    {
        fout << px(x) << "," << py(x) << "," << pz(x) << std::endl;
    }

    fout.close();
    auto bsplineball = curve.BSplineBall();
    filename         = "./result/E_3_Ball" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = t[0]; x <= t[N]; x += 1e-2)
    {
        std::vector<double> p = bsplineball(x);
        fout << p[0] << "," << p[1] << "," << p[2] << std::endl;
    }
    fout.close();

    auto bsplineballproj = curve.BSplineBallProj();
    filename             = "./result/E_3_BallProj" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = t[0]; x <= t[N]; x += 1e-2)
    {
        std::vector<double> p = bsplineballproj(x);
        fout << p[0] << "," << p[1] << "," << p[2] << std::endl;
    }
    fout.close();

    auto ppolyball = curve.PPolyBall();
    filename       = "./result/E_3_PPBall" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = t[0]; x <= t[N]; x += 1e-2)
    {
        std::vector<double> p = ppolyball(x);
        fout << p[0] << "," << p[1] << "," << p[2] << std::endl;
    }
    fout.close();

    auto ppolyballproj = curve.PPolyBallProj();
    filename           = "./result/E_3_PPBallProj" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (double x = t[0]; x <= t[N]; x += 1e-2)
    {
        std::vector<double> p = ppolyballproj(x);
        fout << p[0] << "," << p[1] << "," << p[2] << std::endl;
    }
    fout.close();

    filename = "./result/EXACT" + std::to_string(N) + ".txt";
    fout.open(filename);
    for (int i = 0; i < N + 1; ++i)
    {
        fout << xy[0][i] << "," << xy[1][i] << "," << xy[2][i] << std::endl;
    }
    fout.close();

    Curve<3>      curve_cum(t_cum, xy);
    auto          bspline_cum  = curve_cum.BSpline3d();
    std::string   filename_cum = "./result/E_3_B_cum" + std::to_string(N) + ".txt";
    std::ofstream fout_cum(filename_cum);
    for (double x = t_cum[0]; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << bspline_cum[0](x) << "," << bspline_cum[1](x) << "," << bspline_cum[2](x) << std::endl;
    }

    fout_cum.close();
    auto  ppoly_cum = curve_cum.PPoly3d();
    PPoly px_cum    = ppoly_cum[0];
    PPoly py_cum    = ppoly_cum[1];
    PPoly pz_cum    = ppoly_cum[2];
    filename_cum    = "./result/E_3_PP_cum" + std::to_string(N) + ".txt";
    fout_cum.open(filename_cum);
    for (double x = t_cum[0]; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << px_cum(x) << "," << py_cum(x) << "," << pz_cum(x) << std::endl;
    }
    fout_cum.close();
    auto bsplineball_cum = curve_cum.BSplineBall();
    filename_cum         = "./result/E_3_Ball_cum" + std::to_string(N) + ".txt";
    fout_cum.open(filename_cum);
    for (double x = t_cum[0]; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << bsplineball_cum(x)[0] << "," << bsplineball_cum(x)[1] << "," << bsplineball_cum(x)[2] << std::endl;
    }
    fout_cum.close();
    auto bsplineballproj_cum = curve_cum.BSplineBallProj();
    filename_cum             = "./result/E_3_BallProj_cum" + std::to_string(N) + ".txt";
    fout_cum.open(filename_cum);
    for (double x = t_cum[0]; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << bsplineballproj_cum(x)[0] << "," << bsplineballproj_cum(x)[1] << "," << bsplineballproj_cum(x)[2]
                 << std::endl;
    }
    fout_cum.close();
    auto ppolyball_cum = curve_cum.PPolyBall();
    filename_cum       = "./result/E_3_PPBall_cum" + std::to_string(N) + ".txt";
    fout_cum.open(filename_cum);
    for (double x = t_cum[0]; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << ppolyball_cum(x)[0] << "," << ppolyball_cum(x)[1] << "," << ppolyball_cum(x)[2] << std::endl;
    }
    fout_cum.close();
    auto ppolyballproj_cum = curve_cum.PPolyBallProj();
    filename_cum           = "./result/E_3_PPBallProj_cum" + std::to_string(N) + ".txt";
    fout_cum.open(filename_cum);
    for (double x = t_cum[0]; x <= t_cum[N]; x += 1e-2)
    {
        fout_cum << ppolyballproj_cum(x)[0] << "," << ppolyballproj_cum(x)[1] << "," << ppolyballproj_cum(x)[2]
                 << std::endl;
    }
    fout_cum.close();
    return;
}

int
main(void)
{
    std::cout << "======================E======================" << std::endl;
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
    runB(10);
    runB(40);
    runB(160);
    runC(10);
    runC(40);
    runC(160);
    std::cout << "======================E======================" << std::endl;
}