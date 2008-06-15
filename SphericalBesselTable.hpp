#ifndef __SPHERICALBESSELTABLE_HPP
#define __SPHERICALBESSELTABLE_HPP

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

#include "Defs.hpp"


class SphericalBesselTable
{

    typedef std::vector<Real> RealVector;
    typedef std::vector<RealVector> RealTable2;

    typedef UnsignedInteger Index;


public:

    SphericalBesselTable( const UnsignedInteger nmax,
                          const UnsignedInteger resolution )
	:
	nmax( nmax ),
        delta( M_PI / resolution ),
        sjTable( nmax+1 ),
        syTable( nmax+1 )
    {
	fillTables();
    }

    ~SphericalBesselTable()
    {
        ; // do nothing
    }

    const Real j( const UnsignedInteger n, const Real z ) const
    {
        Real value;

        const RealVector& sjTablen( this->sjTable[n] );
        const Real minz( this->delta * 3 );
        const Real maxz( this->delta * ( sjTablen.size() - 3 ) );

        if( z < minz )
        {
            value = this->_j( n, z );
        }
        else if( z < maxz )
        {
            // table
            value = interp( z, sjTablen, this->delta );
        }
        else
        {
            value = this->_j( n, z );
        }

        return value;
    }

    const Real y( const UnsignedInteger n, const Real z ) const
    {
        Real value;

        const RealVector& syTablen( this->syTable[n] );
        const Real maxz( this->delta * syTablen.size() - 2 );

        if( z < maxz )
        {
            // table
            value = interp( z, syTablen, this->delta );
        }
        else
        {
            value = this->_y( n, z );
        }

        return value;
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

    static const Real maxz( const UnsignedInteger n )
    {
        const Real realn( static_cast<Real>( n ) + 0.5 );

        const Real z( realn * realn + realn + 1 );

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

    void fillTables()
    {
        prepareTables();

        fillzTable();
        fillsjTable();
        fillsyTable();
    }

    void prepareTables()
    {
        for( UnsignedInteger n( 0 ); n <= this->nmax; ++n )
        {
            const Real maxz( this->maxz( n ) );
            const UnsignedInteger 
                maxi( static_cast<UnsignedInteger>( ceil( maxz / delta ) ) );
            this->sjTable[n].resize( maxi );
            this->syTable[n].resize( maxi );
        }
    }

    void fillzTable()
    {
        const UnsignedInteger maxi( this->sjTable[this->nmax].size() );
        this->zTable.reserve( maxi );
        for( UnsignedInteger i( 0 ); i < maxi; ++i )
        {
            const Real z( delta * i );
            this->zTable.push_back( z );
        }
    }

    void fillsjTable()
    {
        RealVector tmpVector( nmax );

        const UnsignedInteger maxi( this->sjTable[this->nmax].size() );

        for( UnsignedInteger i( 0 ); i < maxi; ++i )
        {
            const Real z( delta * i );

            gsl_sf_bessel_jl_steed_array( nmax, z, &tmpVector[0] );

            UnsignedInteger n( nmax );
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
    }

    void fillsyTable()
    {
        RealVector tmpVector( nmax );

        const UnsignedInteger maxi( this->sjTable[this->nmax].size() );
        for( UnsignedInteger i( 1 ); i < maxi; ++i )
        {
            const Real z( delta * i );

            gsl_sf_bessel_yl_array( nmax, z, &tmpVector[0] );

            UnsignedInteger n( nmax );
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
    }

    static const Real interp_simple ( const Real x, const RealVector& table, 
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

    const Real interp( const Real x, const RealVector& table, 
                              const Real delta ) const
    {
       gsl_interp_accel *acc = gsl_interp_accel_alloc();
       const gsl_interp_type *t = gsl_interp_cspline;
       gsl_spline *spline = gsl_spline_alloc(t, table.size());
     
       gsl_spline_init (spline, &this->zTable[0], &table[0], table.size());
     
       const Real y( gsl_spline_eval(spline, x, acc) );
       
       gsl_spline_free (spline);
       gsl_interp_accel_free (acc);

       return y;
    }


private:

    //static const UnsignedInteger n_offset = 2;

    const UnsignedInteger nmax;

    const Real delta;

    //const Real factor;
    //const Real factorp;

    RealVector zTable;

    RealTable2 sjTable;
    RealTable2 syTable;

};



#endif /* __SPHERICALBESSELTABLE_HPP */
