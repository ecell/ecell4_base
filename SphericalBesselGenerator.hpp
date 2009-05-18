#ifndef __SPHERICALBESSELGENERATOR_HPP
#define __SPHERICALBESSELGENERATOR_HPP

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

#include "Defs.hpp"


class SphericalBesselGenerator
{

    typedef UnsignedInteger Index;

public:

    SphericalBesselGenerator()
    {
        ; // do nothing
    }

    ~SphericalBesselGenerator()
    {
        ; // do nothing
    }

    const Real j(const UnsignedInteger n, const Real z) const;

    const Real y(const UnsignedInteger n, const Real z) const;

    static const UnsignedInteger getMinNJ();
    static const UnsignedInteger getMinNY();
    static const UnsignedInteger getMaxNJ();
    static const UnsignedInteger getMaxNY();

    static const SphericalBesselGenerator& instance();

private:

    static const Real _j(const UnsignedInteger n, const Real z)
    {
        return gsl_sf_bessel_jl(n, z);
    }

    static const Real _y(const UnsignedInteger n, const Real z)
    {
        return gsl_sf_bessel_yl(n, z);
    }

    const Real _j_table(const UnsignedInteger n, const Real z) const;
    const Real _y_table(const UnsignedInteger n, const Real z) const;

    static const Real _j_smalln(const UnsignedInteger n, const Real z);
    static const Real _y_smalln(const UnsignedInteger n, const Real z);
     
};




#endif /* __SPHERICALBESSELGENERATOR_HPP */
