#include "BSpline.hpp"
#include "Curve.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
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

int
main(void)
{
    /**
     * @details
     * For the heart curve, we know the equation is
     * \f$ x^2 + (\frac{3}{2}y - \sqrt{|x|})^2 = 3. \f$
     * then we can get
     */
}