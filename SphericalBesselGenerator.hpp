#ifndef __SPHERICALBESSELGENERATOR_HPP
#define __SPHERICALBESSELGENERATOR_HPP

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

#include "Defs.hpp"


class SphericalBesselGenerator
{

    typedef std::vector<Real> RealVector;
    typedef std::vector<RealVector> RealTable2;

    typedef gsl_interp Interpolator;

    typedef std::vector<Interpolator*> InterpolatorVector;

    typedef UnsignedInteger Index;


public:

    SphericalBesselGenerator( const UnsignedInteger maxn,
                              const UnsignedInteger resolution )
	:
	maxn( maxn ),
        delta( M_PI / resolution ),
        interpType( gsl_interp_cspline ),
//        interpType( gsl_interp_akima ),
        interpMinSize( 5 ), // gsl_interp_min_size( ) ),
        sjTable( maxn+1 ),
        syTable( maxn+1 )
    {
	fillTables();
    }

    ~SphericalBesselGenerator()
    {
        ; // do nothing
    }

    const Real j( const UnsignedInteger n, const Real z ) const
    {
        if( n > this->maxn )
        {
            return this->_j( n, z );
        }

        const RealVector& sjTablen( this->sjTable[n] );

        const Real minz( 1.0 + n );
        const Real maxz( this->delta * ( sjTablen.size() - 3 ) );

        if( z < minz )
        {
            return this->_j( n, z );
        }
        else if( z < maxz )
        {
            // table
            return this->_j_table( n, z );
        }
        else
        {
            return this->_j( n, z );
        }
    }

    const Real y( const UnsignedInteger n, const Real z ) const
    {
        if( n > this->maxn )
        {
            return this->_y( n, z );
        }

        const RealVector& syTablen( this->syTable[n] );
        const Real minz( n + 1 );
        const Real maxz( this->delta * ( syTablen.size() - 3 ) );

        if( z < minz )
        {
            return this->_y( n, z );
        }
        else if( z < maxz )
        {
            // table
            return this->_y_table( n, z );
        }
        else
        {
            return this->_y( n, z );
        }
    }

/*
    const Real jp( const int order ) const
    {
	const int index( order-nmin );
	const Real sj( sjArray[index] );
	const Real sjp( ( order * x_r ) * sj - sjArray[index+1] );

	return sjp * factor + sj * factorp;
    }

    const Real yp( const int order ) const
    {
	const int index( order-nmin );
	const Real sy( syArray[index] );
	const Real syp( ( order * x_r ) * sy - syArray[index+1] );

	return syp * factor + sy * factorp;
    }
*/

    /*
    const Real minz( const UnsignedInteger n ) const
    {
        return M_PI;
    }
    */

    const Real maxz( const UnsignedInteger n ) const
    {
        const Real realn( static_cast<Real>( n ) + 0.5 );

        Real z( ( realn * realn + realn + 1 ) * 2e4 );
    
        if( z > 1000. )
        {
            z = std::max( 1000., realn * realn );
        }

        z += this->delta * 10;

        return z;
    }


private:

    static const Real _j( const UnsignedInteger n, const Real z )
    {
        return gsl_sf_bessel_jl( n, z );
    }

    static const Real _y( const UnsignedInteger n, const Real z )
    {
        return gsl_sf_bessel_yl( n, z );
    }

    const Real _j_table( const UnsignedInteger n, const Real z ) const
    {
        return interp( this->sjInterpolatorVector[n],
                       this->zTable,
                       this->sjTable[n],
                       z );
    }

    const Real _y_table( const UnsignedInteger n, const Real z ) const
    {
        return interp( this->syInterpolatorVector[n],
                       this->zTable,
                       this->syTable[n],
                       z );
    }


    void fillTables()
    {
        prepareTables();

        fillzTable();
        fillsjTable();
        fillsyTable();
    }

    void prepareTables()
    {
        for( UnsignedInteger n( 0 ); n <= this->maxn; ++n )
        {
            const Real maxz( this->maxz( n ) );
            const UnsignedInteger 
                maxi( std::max( this->interpMinSize,
                                static_cast<UnsignedInteger>
                                ( ceil( maxz / delta ) ) ) );
            this->sjTable[n].resize( maxi + 1 );
            this->syTable[n].resize( maxi + 1 );
        }
    }

    void fillzTable()
    {
        const UnsignedInteger maxi( this->sjTable[this->maxn].size() - 1 );
        this->zTable.reserve( maxi + 1 );
        for( UnsignedInteger i( 0 ); i <= maxi; ++i )
        {
            const Real z( delta * ( i + 1 ) );
            this->zTable.push_back( z );
        }
    }

    void fillsjTable()
    {
        RealVector tmpVector( maxn + 1 );

        const UnsignedInteger maxi( this->sjTable[this->maxn].size() - 1 );

        for( UnsignedInteger i( 0 ); i <= maxi; ++i )
        {
            const Real z( zTable[i] );

            gsl_sf_bessel_jl_steed_array( maxn, z, &tmpVector[0] );

            UnsignedInteger n( maxn );
            while( this->sjTable[n].size() > i )
            {
                this->sjTable[n][i] = tmpVector[n];

                if( n == 0 )
                {
                    break;
                }

                --n;
            }
        }

        fillInterpolatorVector( this->sjInterpolatorVector, sjTable );
    }

    void fillsyTable()
    {
        RealVector tmpVector( maxn + 1 );

        const UnsignedInteger maxi( this->sjTable[this->maxn].size() - 1 );
        for( UnsignedInteger i( 0 ); i <= maxi; ++i )
        {
            const Real z( zTable[i] );

            gsl_sf_bessel_yl_array( maxn, z, &tmpVector[0] );

            UnsignedInteger n( maxn );
            while( this->syTable[n].size() > i )
            {
                this->syTable[n][i] = tmpVector[n];

                if( n == 0 )
                {
                    break;
                }

                --n;
            }
        }

        fillInterpolatorVector( this->syInterpolatorVector, syTable );
    }

    void fillInterpolatorVector( InterpolatorVector& interpolatorVector,
                                 const RealTable2& table )
    {
        interpolatorVector.resize( this->maxn + 1 );
        for( UnsignedInteger n( 0 ); n <= this->maxn; ++n )
        {
            const RealVector& tablen( table[n] ); 
            gsl_interp* interp = gsl_interp_alloc( this->interpType,
                                                   tablen.size() );

            gsl_interp_init( interp, &this->zTable[0], &tablen[0], 
                             tablen.size() );

            interpolatorVector[n] = interp;

            //printf("n i %d %d\n",n,tablen.size());
        }
    }

    static const Real interp_simple( const Real x, const RealVector& table, 
                                     const Real delta )
    {
        const UnsignedInteger 
            i( static_cast<UnsignedInteger>( floor( x / delta ) ) );

        const Real x0( i * delta );
        const Real x1( ( i + 1 ) * delta );
        
        const Real y0( table[i] );
        const Real y1( table[i+1] );
        
        return y0 + ( x - x0 ) * ( y1 - y0 ) / ( x1 - x0 );
    }

     
    static const Real interp( const Interpolator* interpolator, 
                              const RealVector& xTable,
                              const RealVector& yTable,
                              const Real x )
    {
        //const Real y( gsl_interp_eval( interpolator, &xTable[0],
        //&yTable[0], x, 0 ) );
        //return y;

        return my_cspline_eval( interpolator->state, xTable[1]-xTable[0], 
                                &yTable[0], 
                                yTable.size(), x );
    }

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

    static const double
    my_cspline_eval(const void * vstate,
                    const double dx, const double y_array[], size_t size,
                    double x )
    {
        const cspline_state_t *state = (const cspline_state_t *) vstate;

        //double x_lo, x_hi;
        size_t index;
        
        index = static_cast<size_t>( x / dx ) - 1;
        
        double x_lo = ( index + 1 ) * dx;

        /* evaluate */
        
        const double y_lo = y_array[index];
        const double y_hi = y_array[index + 1];
        const double dy = y_hi - y_lo;
        const double delx = x - x_lo;
        double b_i, c_i, d_i; 
        coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
        return y_lo + delx * (b_i + delx * (c_i + delx * d_i));
    }


private:

    //static const UnsignedInteger n_offset = 2;

    const UnsignedInteger maxn;

    const Real delta;

    const gsl_interp_type* interpType;

    const UnsignedInteger interpMinSize;


    //const Real factor;
    //const Real factorp;

    RealVector zTable;

    InterpolatorVector sjInterpolatorVector;
    InterpolatorVector syInterpolatorVector;
    RealTable2 sjTable;
    RealTable2 syTable;

};




extern SphericalBesselGenerator* sphericalBesselGenerator;

#endif /* __SPHERICALBESSELGENERATOR_HPP */
