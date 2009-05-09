#ifndef __SPHERICALBESSELGENERATOR_HPP
#define __SPHERICALBESSELGENERATOR_HPP

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

#include "Defs.hpp"


class SphericalBesselGenerator
{

    typedef gsl_interp Interpolator;

    typedef std::vector<Interpolator*> InterpolatorVector;

    typedef UnsignedInteger Index;


public:

    SphericalBesselGenerator()
	:
        interpType( gsl_interp_cspline ),
        //interpType( gsl_interp_akima ),
        interpMinSize( 5 ) // gsl_interp_min_size( ) ),
    {
	fillTables();
    }

    ~SphericalBesselGenerator()
    {
        ; // do nothing
    }

    const Real j( const UnsignedInteger n, const Real z ) const;

    const Real y( const UnsignedInteger n, const Real z ) const;

    static const UnsignedInteger getMinNJ();
    static const UnsignedInteger getMinNY();
    static const UnsignedInteger getMaxNJ();
    static const UnsignedInteger getMaxNY();

    static const SphericalBesselGenerator& instance();

private:

    void fillTables();


    /*
    static const sb_table::Table* getSJTable( const UnsignedInteger n );
    static const sb_table::Table* getSYTable( const UnsignedInteger n );
    */

    static const Real _j( const UnsignedInteger n, const Real z )
    {
        return gsl_sf_bessel_jl( n, z );
    }

    static const Real _y( const UnsignedInteger n, const Real z )
    {
        return gsl_sf_bessel_yl( n, z );
    }

    const Real _j_table( const UnsignedInteger n, const Real z ) const;
    const Real _y_table( const UnsignedInteger n, const Real z ) const;

    static const Real _j_smalln( const UnsignedInteger n, const Real z );
    static const Real _y_smalln( const UnsignedInteger n, const Real z );
     
private:

    const gsl_interp_type* interpType;

    const UnsignedInteger interpMinSize;

    InterpolatorVector sjInterpolatorVector;
    InterpolatorVector syInterpolatorVector;

};




#endif /* __SPHERICALBESSELGENERATOR_HPP */
