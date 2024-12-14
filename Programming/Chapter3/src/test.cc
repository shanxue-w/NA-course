#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"
#include "PPoly.hpp"
#include <algorithm>
#include <cmath>
#include <gmpxx.h>

int
main(void)
{
    int N = 100;
    mpf_set_default_prec(2000);
    std::vector<mpf_class> x(N);
    std::vector<mpf_class> y(N);
    for (double i = 0; i < N; i++)
    {
        x[i] = 2.0 * M_PI * i / (N - 1);
        y[i] = sin(x[i].get_d());
    }

    PPInterpolate<4, mpf_class> inter(x, y, 0);
    PPoly<mpf_class>            poly    = inter.getPoly();
    mpf_class                   max_err = mpf_class(0.0);
    for (double i = 0; i < 2 * M_PI; i += 1e-3)
    {
        mpf_class err = abs(inter(i) - sin(i));
        if (err > max_err)
        {
            max_err = err;
        }
    }
    std::cout << "PP-form max error: " << max_err << std::endl;
    BInterpolate<4, mpf_class> inter2(x, y, 0);
    max_err = mpf_class(0.0);
    for (double i = 0; i < 2 * M_PI; i += 1e-3)
    {
        mpf_class err = abs(inter2(i) - std::sin(i));
        if (err > max_err)
        {
            max_err = err;
        }
    }
    std::cout << "BSpline max error: " << max_err << std::endl;
    return 0;
}