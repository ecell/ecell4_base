#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cassert>

#include "compat.h"
#include "SphericalBesselTable.hpp"
#include "SphericalBesselGenerator.hpp"


static inline double hermite_interp(double x, 
                                    double x0, double dx, 
                                    double const* y_array)
{
    const double hinv = 1.0 / dx;

    const size_t i = static_cast<size_t>((x - x0 ) * hinv);
    const size_t index = i * 2;

    const double x_lo = (x - x0) * hinv - i;
    const double x_hi =  1.0 - x_lo;

    const double y_lo = y_array[index];
    const double ydot_lo = y_array[index + 1] * dx;
    const double y_hi = y_array[index + 2];
    const double ydot_hi = y_array[index + 3] * dx;
    
    return x_hi * x_hi * (y_lo + x_lo * (2 * y_lo + ydot_lo)) 
        + x_lo * x_lo * (y_hi + x_hi * (2 * y_hi - ydot_hi));
}

inline static Real interp(Real x_start, Real delta_x,
                          Real const* yTable, Real x)
{
    return hermite_interp(x, x_start, delta_x, yTable);
}

static Real _j(UnsignedInteger n, Real z)
{
    return gsl_sf_bessel_jl(n, z);
}

static Real _y(UnsignedInteger n, Real z)
{
    return gsl_sf_bessel_yl(n, z);
}

SphericalBesselGenerator const& SphericalBesselGenerator::instance()
{
    static const SphericalBesselGenerator sphericalBesselGenerator;
    return sphericalBesselGenerator;
}



UnsignedInteger SphericalBesselGenerator::getMinNJ()
{
    return sb_table::sj_table_min;
}

UnsignedInteger SphericalBesselGenerator::getMinNY()
{
    return sb_table::sy_table_min;
}

UnsignedInteger SphericalBesselGenerator::getMaxNJ()
{
    return sb_table::sj_table_max;
}

UnsignedInteger SphericalBesselGenerator::getMaxNY()
{
    return sb_table::sy_table_max;
}

static sb_table::Table const* getSJTable(UnsignedInteger n)
{
    return sb_table::sj_table[n];
}


static sb_table::Table const* getSYTable(UnsignedInteger n)
{
    return sb_table::sy_table[n];
}

static inline Real _j_table(UnsignedInteger n, Real z)
{
    sb_table::Table const* tablen(getSJTable(n));

    return interp(tablen->x_start, tablen->delta_x, tablen->y, z);
}

static inline Real _y_table(UnsignedInteger n, Real z)
{
    sb_table::Table const* tablen(getSYTable(n));

    return interp(tablen->x_start, tablen->delta_x, tablen->y, z);
}

static inline Real _j_smalln(UnsignedInteger n, Real z)
{
    assert(n <= 3 && n >= 0);

    if(n == 0)
    {
        if(z != 0)
        {
            return std::sin(z) / z;
        }
        else
        {
            return 1.0;
        }
    }

    if(z == 0.0)
    {
        return 0.0;
    }

    Real sin_z;
    Real cos_z;
    sincos(z, &sin_z, &cos_z);

    const Real z_r(1. / z);
        
    if(n == 1)
    {
        return (sin_z * z_r - cos_z) * z_r;
    }
    else if(n == 2)
    {
        const Real _3_zsq(3. * z_r * z_r);
        return (_3_zsq - 1) * sin_z * z_r - _3_zsq * cos_z;
    }
    else //if(n == 3)
    {
        const Real _15_zsq(15. * z_r * z_r);
        return ((_15_zsq - 6.) * sin_z * z_r - 
                (_15_zsq - 1) * cos_z) * z_r;
    }

}

static inline Real _y_smalln(UnsignedInteger n, Real z)
{
    assert(n <= 2 && n >= 0);

    if(n == 0)
    {
        return - std::cos(z) / z;
    }

    Real sin_z;
    Real cos_z;
    sincos(z, &sin_z, &cos_z);

    const Real z_r(1. / z);
        
    if(n == 1)
    {
        return - (cos_z * z_r + sin_z) * z_r;
    }
    else //if(n == 2)
    {
        const Real _3_zsq(3. * z_r * z_r);
        return (1 - _3_zsq) * cos_z * z_r - _3_zsq * sin_z;
    }
}



Real SphericalBesselGenerator::j(UnsignedInteger n, Real z) const
{
    if(n <= 3)
    {
        return _j_smalln(n, z);
    }

    if(n > getMaxNJ())
    {
        return _j(n, z);
    }
    
    const sb_table::Table* table(getSJTable(n));
    assert(table != 0);

    const Real minz(table->x_start + table->delta_x * 3);
    const Real maxz(table->x_start + table->delta_x * (table->N-3));
    
    if(z >= minz && z < maxz)
    {
        return _j_table(n, z);
    }
    else
    {
        return _j(n, z);
    }
}

Real SphericalBesselGenerator::y(const UnsignedInteger n, const Real z) const
{
    if(n <= 2)
    {
        return _y_smalln(n, z);
    }

    if(n > getMaxNY())
    {
        return _y(n, z);
    }
    
    const sb_table::Table* table(getSYTable(n));
    assert(table != 0);
    
    const Real minz(table->x_start + table->delta_x * 3);
    const Real maxz(table->x_start + table->delta_x * (table->N-3));
    
    if(z >= minz && z < maxz)
    {
        return _y_table(n, z);
    }
    else
    {
        return _y(n, z);
    }
}

