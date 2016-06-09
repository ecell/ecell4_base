#include <cmath>
#include "sph_bessel.hpp"

std::pair<std::vector<double>, std::vector<double> > sphj_array(const unsigned int n, const double x) {
    std::vector<double> js(n+1, 0.0), dots(n+1, 0.0);

    if (x == 0) {
        js[0] = 1.0;
        if (n > 0)
            dots[1] = 1.0/3.0;
        return std::make_pair(js, dots);
    }

    js[0] = sin(x)/x;
    dots[0] = cos(x)-js[0]/x;

    if (n == 0)
        return std::make_pair(js, dots);

    js[1] = (js[0]-cos(x))/x;

    for (int k(2); k <= n; ++k)
        js[k] = (2*k-1)*js[k-1]/x - js[k-2];
    for (int k(1); k <= n; ++k)
        dots[k] = js[k-1] - (k+1)*js[k]/x;

    return std::make_pair(js, dots);
}

std::pair<std::vector<double>, std::vector<double> > sphy_array(const unsigned int n, const double x) {
    std::vector<double> ys(n+1, -inf), dots(n+1, inf);

    if (x == 0) {
        return std::make_pair(ys, dots);
    }

    ys[0] = -cos(x)/x;
    dots[0] = (sin(x)-ys[0])/x;

    if (n == 0)
        return std::make_pair(ys, dots);

    ys[1] = (ys[0]-sin(x))/x;

    for (int k(2); k <= n; ++k)
        ys[k] = (2*k-1)*ys[k-1]/x - ys[k-2];
    for (int k(1); k <= n; ++k)
        dots[k] = ys[k-1] - (k+1)*ys[k]/x;

    return std::make_pair(ys, dots);
}
