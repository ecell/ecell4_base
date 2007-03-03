#ifndef __HALFORDERBESSEL_HPP
#define __HALFORDERBESSEL_HPP

#include <cmath>

#include <gsl/gsl_sf_bessel.h>



class HalfOrderBesselGenerator
{

public:

    HalfOrderBesselGenerator( const double x, const int nmin, const int nmax )
	:
	nmin( nmin ),
	x( x ),
	x_r( 1.0 / x ),
	factor( sqrt( ( x + x ) * M_1_PI ) ),
	factorp( sqrt( 1.0 / ( 2.0 * M_PI * x ) ) ),
	//    sjArray( nmax-nmin+2 ),
	//    syArray( nmax-nmin+2 ),
	sjArray( new double[nmax-nmin+2] ),
	syArray( new double[nmax-nmin+2] )
    {
	// nmax+1 order is needed to calculate derivatives.
	fillArrays( nmin, nmax+1 );
    }

    ~HalfOrderBesselGenerator()
    {
	delete[] syArray;
	delete[] sjArray;
    }

    const double j( const int order ) const
    {
	return sjArray[order-nmin] * factor;
    }

    const double y( const int order ) const
    {
	return syArray[order-nmin] * factor;
    }

    const double jp( const int order ) const
    {
	const int index( order-nmin );
	const double sj( sjArray[index] );
	const double sjp( ( order * x_r ) * sj - sjArray[index+1] );

	return sjp * factor + sj * factorp;
    }

    const double yp( const int order ) const
    {
	const int index( order-nmin );
	const double sy( syArray[index] );
	const double syp( ( order * x_r ) * sy - syArray[index+1] );

	return syp * factor + sy * factorp;
    }


    static const double j_0( const double x )
    {
	return std::sin( x ) / x;
    }
    
    static const double y_0( const double x )
    {
	return - std::cos( x ) / x;
    }


private:

    void fillArrays( const int nmin, const int nmax )
    {
	{
	    double jp1( gsl_sf_bessel_jl( nmax+1, x ) );
	    double j( gsl_sf_bessel_jl( nmax, x ) );

	    sjArray[nmax-nmin] = j;

	    for( int n( nmax ); n >= nmin+1; --n )
	    {
		const double jm1( ( n + n + 1.0 ) * x_r * j - jp1 );

		sjArray[n-nmin-1] = jm1;
		jp1 = j;
		j = jm1;
	    }
	}

	{
	    double y( gsl_sf_bessel_yl( nmin+1, x ) );
	    double ym1( gsl_sf_bessel_yl( nmin, x ) );
	    syArray[0] = ym1;
	    syArray[1] = y;


	    for( int n( nmin+1 ); n < nmax; ++n )
	    {
		const double yp1( ( n + n + 1.0 ) * x_r * y - ym1 );

		syArray[n-nmin+1] = yp1;
		ym1 = y;
		y = yp1;
	    }
	}
    }


private:

    const int nmin;
    const double x;
    const double x_r;
    const double factor;
    const double factorp;

    //  DoubleVector sjArray;
    //  DoubleVector syArray;
    double* const sjArray;
    double* const syArray;


};



#endif /* __HALFORDERBESSEL_HPP */
