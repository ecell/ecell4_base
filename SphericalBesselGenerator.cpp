#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "SphericalBesselTable.hpp"

#include "SphericalBesselGenerator.hpp"


/* The modified version of cspline_eval below assumes
   constant interval of samples, and skips bisection search
   at the beginning of the interpolation.  */

/* Taken and adopted from interpolation/cspline.c in GSL */


typedef struct
{
  double * c;
  double * g;
  double * diag;
  double * offdiag;
} cspline_state_t;

static inline
double
my_cspline_eval (const double x,
                 const double x0,
                 const double dx, const double* y_array, 
                 const double* c_array)
{
  const size_t index = trunc((x - x0) / dx);

  const double y_lo = y_array[index];
  const double y_hi = y_array[index + 1];
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];

  const double x_lo = x0 + dx * index;
  const double delx = x - x_lo;

  const double dy = y_hi - y_lo;

  const double b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  const double d = (c_ip1 - c_i) / (3.0 * dx);

  return y_lo + delx * (b + delx * (c_i + delx * d));
}




const SphericalBesselGenerator& SphericalBesselGenerator::instance()
{
    static const SphericalBesselGenerator sphericalBesselGenerator;
    return sphericalBesselGenerator;
}

inline static const Real interp( const gsl_interp* interpolator, 
                                 const Real x_start, const Real delta_x,
                                 const Real* yTable,
                                 const Real x )
{
    return my_cspline_eval(x, x_start, delta_x, yTable,
                           ((cspline_state_t*)interpolator->state)->c); 
}


static void 
fillInterpolatorVector( std::vector<gsl_interp*>& interpolatorVector,
                        const sb_table::Table* const table[],
                        const UnsignedInteger minn,
                        const UnsignedInteger maxn,
                        const gsl_interp_type* interpType )
{
    interpolatorVector.resize( maxn + 1 );
    for( UnsignedInteger n( minn ); n <= maxn; ++n )
    {
        const sb_table::Table* tablen( table[n] ); 
        size_t N(tablen->N);
        gsl_interp* interp = gsl_interp_alloc(interpType, N);
        std::vector<Real> x(N);
        for(size_t i(0); i!= N; ++i)
        {
            x[i] = tablen->x_start + tablen->delta_x * i;
        }
        
        gsl_interp_init( interp, &x[0], tablen->y, tablen->N );
        
        interpolatorVector[n] = interp;
        
        //printf("n i %d %d\n",n,tablen.size());
    }
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

static const sb_table::Table* 
getSJTable( const UnsignedInteger n )
{
    return sb_table::sj_table[n];
}


static const sb_table::Table* 
getSYTable( const UnsignedInteger n )
{
    return sb_table::sy_table[n];
}

inline const Real SphericalBesselGenerator::_j_table( const UnsignedInteger n,
                                                      const Real z ) const
{
    const sb_table::Table* tablen( getSJTable( n ) );

    return interp( this->sjInterpolatorVector[n],
                   tablen->x_start, tablen->delta_x, tablen->y, z );
}

inline const Real SphericalBesselGenerator::_y_table( const UnsignedInteger n, 
                                                      const Real z ) const
{
    const sb_table::Table* tablen( getSYTable( n ) );

    return interp( this->syInterpolatorVector[n],
                   tablen->x_start, tablen->delta_x, tablen->y, z );
}

inline const Real SphericalBesselGenerator::_j_smalln( const UnsignedInteger n,
                                                       const Real z )
{
    assert( n <= 3 && n >= 0 );

    if( n == 0 )
    {
        if( z != 0 )
        {
            return std::sin(z)/z;
        }
        else
        {
            return 1.0;
        }
    }

    if( z == 0.0 )
    {
        return 0.0;
    }

    Real sin_z;
    Real cos_z;
    sincos(z, &sin_z, &cos_z);

    const Real z_r(1. / z);
        
    if(n == 1)
    {
        return ( sin_z * z_r - cos_z ) * z_r;
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

inline const Real SphericalBesselGenerator::_y_smalln( const UnsignedInteger n, 
                                                       const Real z )
{
    assert( n <= 2 && n >= 0 );

    if( n == 0 )
    {
        return - std::cos(z)/z;
    }

    Real sin_z;
    Real cos_z;
    sincos(z, &sin_z, &cos_z);

    const Real z_r(1. / z);
        
    if(n == 1)
    {
        return - ( cos_z * z_r + sin_z ) * z_r;
    }
    else //if(n == 2)
    {
        const Real _3_zsq(3. * z_r * z_r);
        return (1 - _3_zsq) * cos_z * z_r - _3_zsq * sin_z;
    }
}



const Real 
SphericalBesselGenerator::j( const UnsignedInteger n, const Real z ) const
{
    if(n <= 3)
    {
        return _j_smalln(n, z);
    }

    if(n > getMaxNJ())
    {
        return this->_j( n, z );
    }
    
    const sb_table::Table* table( getSJTable( n ) );
    assert( table != 0 );

    const Real minz( table->x_start + table->delta_x * 3 );
    const Real maxz( table->x_start + table->delta_x * (table->N-3) );
    
    if( z >= minz && z < maxz )
    {
        return this->_j_table( n, z );
    }
    else
    {
        return this->_j( n, z );
    }
}

const Real 
SphericalBesselGenerator::y( const UnsignedInteger n, const Real z ) const
{
    if(n <= 2)
    {
        return _y_smalln(n, z);
    }

    if(n > getMaxNY())
    {
        return this->_y( n, z );
    }
    
    const sb_table::Table* table( getSYTable( n ) );
    assert( table != 0 );
    
    const Real minz( table->x_start + table->delta_x * 3 );
    const Real maxz( table->x_start + table->delta_x * (table->N-3) );
    
    if( z >= minz && z < maxz )
    {
        return this->_y_table( n, z );
    }
    else
    {
        return this->_y( n, z );
    }
}

void SphericalBesselGenerator::fillTables()
{
    fillInterpolatorVector( this->sjInterpolatorVector, sb_table::sj_table,
                            getMinNJ(), getMaxNJ(), this->interpType );
    fillInterpolatorVector( this->syInterpolatorVector, sb_table::sy_table,
                            getMinNY(), getMaxNY(), this->interpType );
}

