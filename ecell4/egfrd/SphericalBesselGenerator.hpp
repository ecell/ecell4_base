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

    Real j(UnsignedInteger n, Real z) const;

    Real y(UnsignedInteger n, Real z) const;

    static UnsignedInteger getMinNJ();
    static UnsignedInteger getMinNY();
    static UnsignedInteger getMaxNJ();
    static UnsignedInteger getMaxNY();

    static SphericalBesselGenerator const& instance();
};




#endif /* __SPHERICALBESSELGENERATOR_HPP */
