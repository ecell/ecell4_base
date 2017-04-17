#ifndef ECELL4_SPH_BESSEL_HPP
#define ECELL4_SPH_BESSEL_HPP

#include <utility>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdlib>

typedef std::pair<std::vector<double>, std::vector<double> > values;

const double inf = std::numeric_limits<double>::infinity();

inline double envj(const int n, const double x)
{
    return (0.5*log10(6.28*n)-n*log10(1.36*x/n));
}

inline int nn(int n0, int n1, double f0, double f1, int obj, double x)
{
    int nn;
    double f;
    for (int k(0); k < 20; ++k)
    {
        nn = int(n1 - (n1-n0)/(1-f0/f1));
        f = envj(nn,x)-obj;
        if (abs(nn-n1) < 1)
            break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}

inline int msta1(const double x, const int mp)
{
    const double abs(fabs(x));
    int n0(int(1.1*abs)+1),
        n1(n0+5);
    double f0(envj(n0,abs)-mp),
           f1(envj(n1,abs)-mp);

    return nn(n0, n1, f0, f1, mp, abs);
}

inline int msta2(const double x, const int n, const int mp)
{
    const double abs(fabs(x));
    const double hmp(0.5*mp);
    const double ejn(envj(n,abs));

    double obj;
    int n0;
    if (ejn <= hmp)
    {
        obj = mp;
        n0 = int(1.1*abs)+1;
    }
    else
    {
        obj = hmp + ejn;
        n0 = n;
    }

    double f0(envj(n0,abs)-obj);
    int n1(n0+5);
    double f1(envj(n1,abs)-obj);

    return nn(n0, n1, f0, f1, obj, abs)+10;
}

inline values sphj_array(const int n, const double x)
{
    std::vector<double> js(n+1, 0.0), dots(n+1, 0.0);

    if (x == 0)
    {
        js[0] = 1.0;
        if (n > 0)
            dots[1] = 1.0/3.0;
        return std::make_pair(js, dots);
    }

    js[0] = sin(x)/x;
    dots[0] = (cos(x)-js[0])/x;

    if (n == 0)
        return std::make_pair(js, dots);

    js[1] = (js[0]-cos(x))/x;

    int maxn(n);
    if (n >= 2)
    {
        const double j0(js[0]), j1(js[1]);
        int m(msta1(x, 200));
        if (m < maxn)
            maxn = m;
        else
            m = msta2(x, maxn, 15);

        if (maxn != n)
            throw "sphj_array precision error";

        double f(0.0), f0(0.0), f1(1.0e0-100);
        for (int k(m); k >= 0; --k)
        {
            f = (2*k+3)*f1/x - f0;
            if (k <= maxn)
                js[k] = f;
            f0 = f1;
            f1 = f;
        }

        double c;
        if (fabs(j0) > fabs(j1))
            c = j0/f;
        else
            c = j1/f0;

        for (int k(0); k <= maxn; ++k)
            js[k] = c * js[k];
    }

    for (int k(1); k <= maxn; ++k)
        dots[k] = js[k-1] - (k+1)*js[k]/x;

    return std::make_pair(js, dots);
}

inline values sphy_array(const int n, const double x)
{
    std::vector<double> ys(n+1, -inf), dots(n+1, inf);

    if (x == 0)
        return std::make_pair(ys, dots);

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

#endif /* ECELL4_SPH_BESSEL_HPP */
