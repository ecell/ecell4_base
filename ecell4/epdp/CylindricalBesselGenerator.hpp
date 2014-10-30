#ifndef __CYLINDRICALBESSELGENERATOR_HPP
#define __CYLINDRICALBESSELGENERATOR_HPP

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include "Defs.hpp"


class CylindricalBesselGenerator
{

    typedef UnsignedInteger Index;

public:

    CylindricalBesselGenerator()
    {
        ; // do nothing
    }

    ~CylindricalBesselGenerator()
    {
        ; // do nothing
    }

    Real J(UnsignedInteger n, Real z) const;

    Real Y(UnsignedInteger n, Real z) const;

    static UnsignedInteger getMinNJ();
    static UnsignedInteger getMinNY();
    static UnsignedInteger getMaxNJ();
    static UnsignedInteger getMaxNY();

    static CylindricalBesselGenerator const& instance();
};




#endif /* __CYLINDRICALBESSELGENERATOR_HPP */
