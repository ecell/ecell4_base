#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cassert>

#include "compat.h"
#include "CylindricalBesselTable.hpp"
#include "CylindricalBesselGenerator.hpp"


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

static Real _J(UnsignedInteger n, Real z)
{
    return gsl_sf_bessel_Jn(n, z);
}

static Real _Y(UnsignedInteger n, Real z)
{
    return gsl_sf_bessel_Yn(n, z);
}

CylindricalBesselGenerator const& CylindricalBesselGenerator::instance()
{
    static const CylindricalBesselGenerator cylindricalBesselGenerator;
    return cylindricalBesselGenerator;
}



UnsignedInteger CylindricalBesselGenerator::getMinNJ()
{
    return cb_table::cj_table_min;
}

UnsignedInteger CylindricalBesselGenerator::getMinNY()
{
    return cb_table::cy_table_min;
}

UnsignedInteger CylindricalBesselGenerator::getMaxNJ()
{
    return cb_table::cj_table_max;
}

UnsignedInteger CylindricalBesselGenerator::getMaxNY()
{
    return cb_table::cy_table_max;
}

static cb_table::Table const* getCJTable(UnsignedInteger n)
{
    return cb_table::cj_table[n];
}


static cb_table::Table const* getCYTable(UnsignedInteger n)
{
    return cb_table::cy_table[n];
}

static inline Real _J_table(UnsignedInteger n, Real z)
{
    cb_table::Table const* tablen(getCJTable(n));

    return interp(tablen->x_start, tablen->delta_x, tablen->y, z);
}

static inline Real _Y_table(UnsignedInteger n, Real z)
{
    cb_table::Table const* tablen(getCYTable(n));

    return interp(tablen->x_start, tablen->delta_x, tablen->y, z);
}



Real CylindricalBesselGenerator::J(UnsignedInteger n, Real z) const
{
    if(n > getMaxNJ())
    {
        return _J(n, z);
    }
    
    const cb_table::Table* table(getCJTable(n));
    assert(table != 0);

    const Real minz(table->x_start + table->delta_x * 3);
    const Real maxz(table->x_start + table->delta_x * (table->N-3));
    
    if(z >= minz && z < maxz)
    {
        return _J_table(n, z);
    }
    else
    {
        return _J(n, z);
    }
}

Real CylindricalBesselGenerator::Y(const UnsignedInteger n, const Real z) const
{
    if(n > getMaxNY())
    {
        return _Y(n, z);
    }
    
    const cb_table::Table* table(getCYTable(n));
    assert(table != 0);
    
    const Real minz(table->x_start + table->delta_x * 3);
    const Real maxz(table->x_start + table->delta_x * (table->N-3));
    
    if(z >= minz && z < maxz)
    {
        return _Y_table(n, z);
    }
    else
    {
        return _Y(n, z);
    }
}

