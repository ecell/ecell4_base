#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cassert>

#include "SphericalBesselTable.hpp"

#include "SphericalBesselGenerator.hpp"


static inline double hermite_interp(const double x, 
                                    const double x0, const double dx, 
                                    const double* y_array)
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

inline static const Real interp(const Real x_start, const Real delta_x,
                                const Real* yTable, const Real x)
{
    return hermite_interp(x, x_start, delta_x, yTable);
}

const SphericalBesselGenerator& SphericalBesselGenerator::instance()
{
    static const SphericalBesselGenerator sphericalBesselGenerator;
    return sphericalBesselGenerator;
}



const UnsignedInteger
SphericalBesselGenerator::getMinNJ()
{
    return sb_table::sj_table_min;
}

const UnsignedInteger
SphericalBesselGenerator::getMinNY()
{
    return sb_table::sy_table_min;
}

const UnsignedInteger
SphericalBesselGenerator::getMaxNJ()
{
    return sb_table::sj_table_max;
}

const UnsignedInteger
SphericalBesselGenerator::getMaxNY()
{
    return sb_table::sy_table_max;
}

static const sb_table::Table* getSJTable(const UnsignedInteger n)
{
    return sb_table::sj_table[n];
}


static const sb_table::Table* getSYTable(const UnsignedInteger n)
{
    return sb_table::sy_table[n];
}

inline const Real SphericalBesselGenerator::_j_table(const UnsignedInteger n,
                                                     const Real z) const
{
    const sb_table::Table* tablen(getSJTable(n));

    return interp(tablen->x_start, tablen->delta_x, tablen->y, z);
}

inline const Real SphericalBesselGenerator::_y_table(const UnsignedInteger n, 
                                                     const Real z) const
{
    const sb_table::Table* tablen(getSYTable(n));

    return interp(tablen->x_start, tablen->delta_x, tablen->y, z);
}

inline const Real SphericalBesselGenerator::_j_smalln(const UnsignedInteger n,
                                                      const Real z)
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

inline const Real SphericalBesselGenerator::_y_smalln(const UnsignedInteger n, 
                                                      const Real z)
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



const Real 
SphericalBesselGenerator::j(const UnsignedInteger n, const Real z) const
{
    if(n <= 3)
    {
        return _j_smalln(n, z);
    }

    if(n > getMaxNJ())
    {
        return this->_j(n, z);
    }
    
    const sb_table::Table* table(getSJTable(n));
    assert(table != 0);

    const Real minz(table->x_start + table->delta_x * 3);
    const Real maxz(table->x_start + table->delta_x * (table->N-3));
    
    if(z >= minz && z < maxz)
    {
        return this->_j_table(n, z);
    }
    else
    {
        return this->_j(n, z);
    }
}

const Real 
SphericalBesselGenerator::y(const UnsignedInteger n, const Real z) const
{
    if(n <= 2)
    {
        return _y_smalln(n, z);
    }

    if(n > getMaxNY())
    {
        return this->_y(n, z);
    }
    
    const sb_table::Table* table(getSYTable(n));
    assert(table != 0);
    
    const Real minz(table->x_start + table->delta_x * 3);
    const Real maxz(table->x_start + table->delta_x * (table->N-3));
    
    if(z >= minz && z < maxz)
    {
        return this->_y_table(n, z);
    }
    else
    {
        return this->_y(n, z);
    }
}

