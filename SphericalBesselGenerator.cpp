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

static inline void
coeff_calc (const double c_array[], double dy, double dx, size_t index,  
            double * b, double * c, double * d)
{
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  *c = c_i;
  *d = (c_ip1 - c_i) / (3.0 * dx);
}


static inline
int
my_cspline_eval (const void * vstate,
                 const double x_array[], const double y_array[], size_t size,
                 double x,
                 gsl_interp_accel * a,
                 double *y)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  const double x0 = x_array[0];
  const double dx = x_array[1] - x0;

  const size_t index = trunc((x - x0) / dx);

  //const double x_lo = x_array[index];
  const double x_lo = x0 + dx * index;

  const double y_lo = y_array[index];
  const double y_hi = y_array[index + 1];
  const double dy = y_hi - y_lo;
  const double delx = x - x_lo;
  double b_i, c_i, d_i; 
  coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
  *y = y_lo + delx * (b_i + delx * (c_i + delx * d_i));
  return GSL_SUCCESS;
}




const SphericalBesselGenerator& SphericalBesselGenerator::instance()
{
    static const SphericalBesselGenerator sphericalBesselGenerator;
    return sphericalBesselGenerator;
}

inline static const Real interp( const gsl_interp* interpolator, 
                                 const Real* xTable, const Real* yTable,
                                 const Real x )
{
//    const Real y( gsl_interp_eval( interpolator, xTable, yTable, x, 0 ) );
    
    Real y;
    my_cspline_eval( interpolator->state, xTable, yTable,
                     interpolator->size, x, 0, &y );

    return y;
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
        gsl_interp* interp = gsl_interp_alloc( interpType,
                                               tablen->N );
        
        gsl_interp_init( interp, tablen->x, tablen->y, tablen->N );
        
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
                   tablen->x, tablen->y, z );
}

inline const Real SphericalBesselGenerator::_y_table( const UnsignedInteger n, 
                                                      const Real z ) const
{
    const sb_table::Table* tablen( getSYTable( n ) );

    return interp( this->syInterpolatorVector[n],
                   tablen->x, tablen->y, z );
}


const Real 
SphericalBesselGenerator::j( const UnsignedInteger n, const Real z ) const
{
    if( n > getMaxNJ() || n < getMinNJ() )
    {
        return this->_j( n, z );
    }
    
    const sb_table::Table* table( getSJTable( n ) );
    assert( table != 0 );

    const Real minz( table->x[3] );
    //const Real minz( 1.0 * n );
    const Real maxz( table->x[table->N - 3] );
    
    if( z > minz && z < maxz )
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
    if( n > getMaxNY() || n < getMinNY() )
    {
        return this->_y( n, z );
    }
    
    const sb_table::Table* table( getSYTable( n ) );
    assert( table != 0 );
    
    const Real minz( table->x[3] );
    //const Real minz( 1.0 * n );
    const Real maxz( table->x[table->N - 3] );
    
    if( z > minz && z < maxz )
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

