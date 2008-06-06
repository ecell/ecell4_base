//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>

#include "factorial.hpp"
#include "funcSum.hpp"
#include "findRoot.hpp"

#include "FirstPassagePairGreensFunction.hpp"



FirstPassagePairGreensFunction::
FirstPassagePairGreensFunction( const Real D, 
				const Real kf, 
				const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
//    alphaTable( this->MAX_ORDER+1 ),
    h( kf / ( 4.0 * M_PI * Sigma * Sigma * D ) ),
    hsigma_p_1( 1.0 + h * Sigma ),
    a( INFINITY )
//    alpha0_threshold( 0.0 )
{
    ; // do nothing
}

FirstPassagePairGreensFunction::~FirstPassagePairGreensFunction()
{
    ; // do nothing
}

void FirstPassagePairGreensFunction::seta( const Real a )
{
    const Real sigma( this->getSigma() );

    THROW_UNLESS( std::invalid_argument, a >= sigma );

    this->a = a;

    clearAlphaTable();
}

//
// Alpha-related methods
//

void FirstPassagePairGreensFunction::clearAlphaTable() const
{
    std::for_each( this->alphaTable.begin(), this->alphaTable.end(),
		   boost::mem_fn( &RealVector::clear ) );
    this->alphaOffsetTable[0] = 0;
    std::fill( this->alphaOffsetTable.begin()+1, this->alphaOffsetTable.end(),
	       -1 );

}


const Real 
FirstPassagePairGreensFunction::f_alpha0( const Real alpha ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real alpha_a_m_sigma( alpha * ( a - sigma ) );
    const Real hsigma_p_1( this->hsigma_p_1 );

    Real sin_alpha_a_m_sigma;
    Real cos_alpha_a_m_sigma;
    sincos( alpha_a_m_sigma, &sin_alpha_a_m_sigma, &cos_alpha_a_m_sigma );

    const Real term1( alpha * sigma * cos_alpha_a_m_sigma );
    const Real term2( hsigma_p_1 * sin_alpha_a_m_sigma );

    const Real result( term1 + term2 );

    return result;
}


const Real 
FirstPassagePairGreensFunction::f_alpha0_aux( const Real alpha ) const
{
    const Real a( this->geta() );
    const Real sigma( this->getSigma() );

    const Real term1( ( a - sigma ) * alpha );

    const Real angle( this->hsigma_p_1 / ( sigma * alpha ) );
    const Real term2( std::atan( angle ) );

    const Real result( term1 - term2 );

    return result;
}




const Real 
FirstPassagePairGreensFunction::
f_alpha0_aux_F( const Real alpha, const f_alpha0_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );

    return gf->f_alpha0_aux( alpha ) - value;
//    return gf->f_alpha0( alpha );
}


const Real 
FirstPassagePairGreensFunction::alpha0_i( const Integer i ) const
{
    THROW_UNLESS( std::out_of_range, i >= 0 );

    const Real a( this->geta() );
    const Real sigma( this->getSigma() );


    const Real target( i * M_PI + M_PI_2 );
    f_alpha0_aux_params params = { this, target };


    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &f_alpha0_aux_F ),
	    &params 
	};


    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real interval( M_PI / ( a - sigma ) );
    Real low( i * interval + std::numeric_limits<Real>::epsilon() );
    Real high( (i+1) * interval );

    //assert( GSL_FN_EVAL( &F, low ) * GSL_FN_EVAL( &F, high ) < 0.0 );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int j( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

        low = gsl_root_fsolver_x_lower( solver );
        high = gsl_root_fsolver_x_upper( solver );
	const int status( gsl_root_test_interval( low, high, 0.0, 1e-15 ) );

	if( status == GSL_CONTINUE )
	{
	    if( j >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha0_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++j;
    }

    const Real alpha( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
  
//    printf("alphavalue %g, alpha %g\n",GSL_FN_EVAL( &F, alpha ), alpha);

    return alpha;
}


void
FirstPassagePairGreensFunction::updateAlphaTable0( const Real t ) const
{
    RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
    alphaTable_0.clear();
    alphaTable_0.reserve( MAX_ALPHA_SEQ );

    const Real alpha0_0( this->alpha0_i( 0 ) );
    alphaTable_0.push_back( alpha0_0 );

    const Real Dt( this->getD() * t );

//    const Real alpha_cutoff( sqrt( ( - log( TOLERANCE * 1e-2 ) / Dt )
//				   + alpha0_0 * alpha0_0 ) );
    const Real alpha_cutoff( sqrt( ( - log( TOLERANCE * 1e-3 ) / Dt ) ) );

    unsigned int i( 1 );
    while( true )
    {
	const Real alpha0_i( this->alpha0_i( i ) );
	alphaTable_0.push_back( alpha0_i );

	if( alpha0_i > alpha_cutoff && i >= 10 ) // make at least 10 terms
	{
            //printf("\nalpha n %d\n",i );
	    break;
	}

	++i;

	if( i >= MAX_ALPHA_SEQ )
	{
	    break;
	}
    }
}

const Real FirstPassagePairGreensFunction::f_alpha( const Real alpha,
						    const Integer n ) const
{
    const Real a( this->geta() );
    const Real sigma( getSigma() );
    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( getSigma() * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );

    const Real hSigma_m_n( hSigma - realn );

    // GSL
    const Real jas1( gsl_sf_bessel_jl( n, sigmaAlpha ) );
    const Real yas1( gsl_sf_bessel_yl( n, sigmaAlpha ) );
    const Real jas2( gsl_sf_bessel_jl( n+1, sigmaAlpha ) );
    const Real yas2( gsl_sf_bessel_yl( n+1, sigmaAlpha ) );
    const Real jaa( gsl_sf_bessel_jl( n, aAlpha ) );
    const Real yaa( gsl_sf_bessel_yl( n, aAlpha ) );


//    printf("b %g %g %g %g %g %g\n", jas1, jas2, yas1, yas2, yaa, jaa );
    const Real term1( ( hSigma_m_n * jas1 + sigmaAlpha * jas2 ) * yaa );
    const Real term2( ( hSigma_m_n * yas1 + sigmaAlpha * yas2 ) * jaa );

//    printf("s %g %g %g %g\n", hSigma_m_n * jas1 * yaa, sigmaAlpha * jas2 * yaa,
//	   hSigma_m_n * yas1 * jaa, sigmaAlpha * yas2 * jaa);

//    printf("t %g %g %g %g\n", alpha, term1, term2, term1-term2 );// cos(f_alpha_aux( alpha,n )) );
    const Real factor( 2.0 * alpha * sqrt( a * sigma ) * M_1_PI );

    const Real result( ( term1 - term2 ) * factor );
    
    return result;
}

inline const Real G( const unsigned int n, const unsigned int k )
{
    //    std::cerr << n << ' ' << k << std::endl;
    //    return gsl_sf_fact( n + k ) / ( gsl_sf_fact( k ) * 
    //                                    gsl_sf_fact( n - k ) );    
    return factorial( n + k ) * ( factorial_r( k ) * factorial_r( n - k ) );
}


const Real FirstPassagePairGreensFunction::P( const Integer n,
                                              const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 );
    Integer term1( 1 );

    const Real x2sq_r( 1.0 / gsl_pow_2( x + x ) );
    const unsigned int maxm( n / 2 );
    for( unsigned int m( 0 ); m <= maxm; ++m )
    {
	const Real value( term1 * sx2 * G( n, 2 * m ) );
	result += value;

	term1 = - term1;
	sx2 *= x2sq_r;
    }

    return result;
}

const boost::tuple<Real,Real>
FirstPassagePairGreensFunction::P2( const Integer n, const Real x )
{
    Real result( 0.0 );
    Real resultp( 0.0 );

    Real sx2( 1.0 );
    Integer term1( 1 );

    const Real x2sq_r( 1.0 / gsl_pow_2( x + x ) );
    const unsigned int np1( n + 1 );
    const unsigned int maxm( n / 2 );
    for( unsigned int m( 0 ); m <= maxm; ++m )
    {
	const Real sx2p( term1 * sx2 );
	const unsigned int m2( 2 * m );
	const Real value( sx2p * G( n, m2 ) );
	result += value;

	const Real valuep( sx2p * G( np1, m2 ) );
	resultp += valuep;

	term1 = - term1;
	sx2 *= x2sq_r;
    }

    if( n % 2 )
    {
	resultp += term1 * sx2 * G( np1, np1 );
    }


    return boost::make_tuple( result, resultp );
}


const Real FirstPassagePairGreensFunction::Q( const Integer n,
					      const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 / ( x + x ) );
    Integer term1( 1 );

    const Real x2sq( sx2 * sx2 );
    const unsigned int maxm( (n+1)/2 ); // sum_(0)^((n-1)/2)
    for( unsigned int m( 0 ); m < maxm; ++m )
    {
	const Real value( term1 * sx2 * G( n, 2 * m + 1 ) );
	result += value;

	term1 = - term1;  // (-1)^m
	sx2 *= x2sq;
    }

    return result;
}

const boost::tuple<Real,Real>
FirstPassagePairGreensFunction::Q2( const Integer n, const Real x )
{
    Real result( 0.0 );
    Real resultp( 0.0 );

    Real sx2( 1.0 / ( x + x ) );
    Integer term1( 1 );  // (-1)^m

    const Real x2sq( sx2 * sx2 );
    const unsigned int np1( n + 1 );
    const unsigned int maxm( (n+1)/2 ); // sum_(0)^((n-1)/2)
    for( unsigned int m( 0 ); m < maxm; ++m )
    {
	const Real sx2p( term1 * sx2 );
	const unsigned int m2p1( 2 * m + 1 );
	const Real value( sx2p * G( n, m2p1 ) );
	result += value;

	const Real valuep( sx2p * G( np1, m2p1 ) );
	resultp += valuep;

	term1 = - term1; // (-1)^m
	sx2 *= x2sq;
    } 


    if( !( n % 2 ) )
    {
	resultp += term1 * sx2 * G( np1, np1 );
    }


    return boost::make_tuple( result, resultp );
}


const Real 
FirstPassagePairGreensFunction::f_alpha_aux( const Real alpha, 
					     const Integer n ) const
{
    if( alpha == 0.0 )
    {
	return -1.0;
    }

    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );

    const Real n_m_hSigma( n - h * sigma );

    /*(a - s) u - 
      ArcTan[( P[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) -
               Q[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]) )/
             ( Q[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) + 
               P[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]) )]
    */

    const Real Pa( P( n, aAlpha ) );
    const Real Qa( Q( n, aAlpha ) );

    Real Ps;
    Real Psp;
    boost::tie( Ps, Psp ) = P2( n, sigmaAlpha );

    Real Qs;
    Real Qsp;
    boost::tie( Qs, Qsp ) = Q2( n, sigmaAlpha );

    const Real n_m_hSigmaPs( n_m_hSigma * Ps );
    const Real n_m_hSigmaQs( n_m_hSigma * Qs );
    const Real sigmaAlphaPsp( sigmaAlpha * Psp );
    const Real sigmaAlphaQsp( sigmaAlpha * Qsp );

    const Real Qa_Pa( Qa / Pa );

    const Real A( sigmaAlphaQsp - n_m_hSigmaPs );
    const Real B( sigmaAlphaPsp + n_m_hSigmaQs );

    // this form, dividing all terms by Pa, prevents overflow.
    const Real angle( ( A - Qa_Pa * B ) / ( Qa_Pa * A + B ) );

    const Real term1( ( a - sigma ) * alpha );
    const Real term2( std::atan( angle ) );

    const Real result( term1 - term2 );

    /*
    if( ! finite( result ) )// debug
    {
        printf("alpha %g n %d\n",alpha,n );
        printf("t1 %g t2 %g t3 %g t4 %g\n",t1,t2,t3,t4);
        printf("Pa %g Qa %g Ps %g Qs %g Psp %g Qsp %g\n",Pa,Qa,Ps,Qs,Psp,Qsp);
        printf("aux %g %g %g %g\n",angle,term1, term2, result );
    }
    */

    return result;
}


const Real 
FirstPassagePairGreensFunction::
f_alpha_aux_F( const Real alpha,
	       const f_alpha_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Integer n( params->n );
    const Real value( params->value );


    return gf->f_alpha_aux( alpha, n ) - value;
//    return gf->f_alpha( alpha, n );
}


const Real 
FirstPassagePairGreensFunction::alpha_i( const Integer i, const Integer n, 
					 gsl_root_fsolver* const solver ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    const Real target( M_PI * i + M_PI_2 );

    const Real factor( 1.0 / ( a - sigma ) );
    Real low( ( target - M_PI_2 ) * factor );
    Real high( ( target + M_PI_2 ) * factor );

    f_alpha_aux_params params = { this, n, target };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>
	    ( &FirstPassagePairGreensFunction::f_alpha_aux_F ),
	    &params 
	};

    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );
    unsigned int k( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	
	low = gsl_root_fsolver_x_lower( solver );
	high = gsl_root_fsolver_x_upper( solver );
	const int status( gsl_root_test_interval( low, high, 1e-6, 1e-15 ) );
	
	if( status == GSL_CONTINUE )
	{
	    if( k >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}
	
	++k;
    }
    
    const Real alpha( gsl_root_fsolver_root( solver ) );

    return alpha;
}


const unsigned int
FirstPassagePairGreensFunction::alphaOffset( const unsigned int n ) const
{
    if( this->alphaOffsetTable[n] >= 0 )
    {
	return this->alphaOffsetTable[n];
    }

    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    assert( this->alphaOffsetTable.size() >= n );
    unsigned int offset( this->alphaOffsetTable[n-1] );

    const Real factor( 1.0 / ( a - sigma ) );

    Real target( offset * M_PI + M_PI_2 );
    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real alphaMid( target * factor );
    const Real alphaHalfRange( M_PI_2 * factor );
    Real low( alphaMid - alphaHalfRange * (1.0 - 1e-3) ); // avoid zero.
    Real high( alphaMid + alphaHalfRange );


    // Here we find the interval where the first positive root is in.
    // We find the first pair of alpha
    // ( Pi * offset + Pi/2 ) +- Pi/2 / ( a - sigma )
    // where the values of f_alpha() straddle.
    // The assumption is the interval between roots is not much
    // smaller than Pi / ( a - sigma ).


    Real lowvalue( f_alpha(low,n) );
    Real highvalue( f_alpha(high,n) );

    while( true ) // this can be much faster if better initial guess is given.
    {

	if( lowvalue * highvalue < 0 ) // low and high straddle?
	{
	    break;
	}

	//printf("lh: %d %d %g %g %g %g\n", 
        //n, offset, low, high, f_alpha( low, n ), f_alpha( high, n ) );
	++offset;
	target = M_PI * offset + M_PI_2;
	low = ( target - M_PI_2 ) * factor;
	high = ( target + M_PI_2 ) * factor;

	lowvalue = highvalue;
	highvalue = f_alpha( high, n );
    }

    this->alphaOffsetTable[n] = offset;

    return offset;
}


void
FirstPassagePairGreensFunction::updateAlphaTable( const unsigned int n,
						  const Real t ) const
{
    THROW_UNLESS( std::range_error, n >= 0 && n <= this->MAX_ORDER );

    if( n == 0 )
    {
	this->updateAlphaTable0( t );
	return;
    }

    const unsigned int offset( alphaOffset( n ) );
//    const unsigned int offset( 0 );

    RealVector& alphaTable_n( this->getAlphaTable( n ) );
    alphaTable_n.clear();
    alphaTable_n.reserve( MAX_ALPHA_SEQ );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real alphan_0( alpha_i( offset, n, solver ) );
    const Real alphan_0_sq( alphan_0 * alphan_0 );

    alphaTable_n.push_back( alphan_0 );

    const Real Dt( this->getD() * t );

    const Real threshold( this->TOLERANCE * 1e-2 * 
			  alphan_0_sq * exp( - Dt * alphan_0_sq ) );
   
    const unsigned int end( offset + MAX_ALPHA_SEQ );
    unsigned int i( offset + 1 );
    while( true )
    {
	const Real alpha_i( this->alpha_i( i, n, solver ) );

//	printf("alpha %d %d %g %g %g\n", n, i, alpha_i, f_alpha(alpha_i+1,n),
//	       f_alpha(alpha_i-1,n));

	alphaTable_n.push_back( alpha_i );

	// cutoff
	const Real alpha_i_sq( alpha_i * alpha_i );
	if( alpha_i_sq * exp( - Dt * alpha_i_sq )  < threshold )
	{
	    break;
	}


	++i;

	if( i >= end )
	{
	    std::cerr << "alphaTable (" << n << 
		"): didn't converge. t = " << t  << ", " 
		      << dump() << std::endl;
	    break;
	}
    }

    gsl_root_fsolver_free( solver );
}




const Real 
FirstPassagePairGreensFunction::p_0_i( const Real alpha, 
				       const Real r,
				       const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    Real num1;
    {
	const Real angle_r( alpha * ( r - sigma ) );
	Real sin_r;
	Real cos_r;
	sincos( angle_r, &sin_r, &cos_r );
	num1 = alpha * sigma * cos_r + hsigma_p_1 * sin_r ;
    }

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( 2 * M_PI * r * r0 * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_survival_i( const Real alpha,
					      const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real angle_a( alpha * ( a - sigma ) );
    const Real cos_a( cos( angle_a ) );

    const Real num1( h * sigmasq * hsigma_p_1 
                     - a * ( hsigma_p_1 * hsigma_p_1
                             + sigmasq * alphasq ) * cos_a );

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( r0 * hsigma_p_1 * alpha * 
                    ( - hsigma_p_1 *
                      ( a + a * h * sigma - h * sigmasq ) 
                      + ( sigma - a ) * sigmasq * alphasq ) );

    const Real result( - 2.0 * num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::dp_survival_i( const Real alpha,
					       const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real angle_a( alpha * ( a - sigma ) );
    const Real cos_a( cos( angle_a ) );

    const Real num1( alpha * ( h * sigmasq * hsigma_p_1 
                               - a * ( hsigma_p_1 * hsigma_p_1 
                                       + sigmasq * alphasq ) *
                               cos_a ) );

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( r0 * hsigma_p_1 * 
                    (- hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) )
                    + ( sigma - a ) * sigmasq * alphasq );

    const Real result( 2.0 * getD() * num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::leavea_i( const Real alpha,
					  const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real D( getD() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );
    const Real hsq( h * h );

    const Real f1( 1 + 2 * h * sigma + hsq * sigmasq + sigmasq * alphasq );

    const Real angle_a( alpha * ( a - sigma ) );
    const Real cos_a( cos( angle_a ) );

    const Real num1( alpha * f1 * cos_a );
    const Real num2( num_r0( alpha, r0 ) );
    
    const Real den( 2 * a * M_PI * r0 * hsigma_p_1 *
                    ( a * f1
                      - sigmasq * ( h + hsq * sigma + sigma * alphasq ) ) );

    const Real result( D * num1 * num2 / den );
    //printf("leavea_i %g %g %g %g\n", result, num1, num2, den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::leaves_i( const Real alpha,
					  const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real D( getD() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real num( h * alpha * num_r0( alpha, r0 ) );
		      
    const Real den( 2 * M_PI * r0 *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

//    printf("leaves_i %g %g %g\n", h * alpha * num_r0( alpha, r0 ),
//          ( a - sigma ) * sigmasq * alphasq,
//           hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) );
    const Real result( - D * num / den );
	
    return result;
}


const Real 
FirstPassagePairGreensFunction::p_leavea_i( const Real alpha,
					    const Real r0,
                                            const Real pleave_factor ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real hsigma_p_1( this->hsigma_p_1 );
    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real angle_a( alpha * ( a - sigma ) );
    const Real cos_a( cos( angle_a ) );

    const Real num1( ( hsigma_p_1 * hsigma_p_1 + sigmasq * alphasq ) * cos_a );

    const Real result( - 2.0 * a * num1 * pleave_factor / hsigma_p_1 );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_leaves_i( const Real alpha,
					    const Real r0,
                                            const Real pleave_factor ) const
{
    const Real sigma( getSigma() );
    const Real h( geth() );
 
    const Real num( h * sigma * sigma );
		      
    const Real result( 2.0 * num * pleave_factor );
	
    return result;
}

const Real 
FirstPassagePairGreensFunction::p_survival_den( const Real alpha,
                                                const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );
    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real den( r0 * alpha *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );
    
    return den;
}



const Real FirstPassagePairGreensFunction::num_r0( const Real alpha,
                                                   const Real r0 ) const
{
    const Real sigma( getSigma() );
    const Real angle_r0( alpha * ( r0 - sigma ) );

    Real sin_r0;
    Real cos_r0;
    sincos( angle_r0, &sin_r0, &cos_r0 );

    const Real hsigma_p_1( this->hsigma_p_1 );
    const Real result( alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0 );

    return result;
}


const Real FirstPassagePairGreensFunction::pleaveFactor( const Real alpha,
                                                         const Real r0 ) const
{
    return num_r0( alpha, r0 ) / p_survival_den( alpha, r0 );
}


const Real
FirstPassagePairGreensFunction::p_int_r_i( const Real r,
					   const Real alpha,
					   const Real r0,
					   const Real num_r0 ) const
{
    const Real sigma( getSigma() );

    const Real angle_r( alpha * ( r - sigma ) );
    Real sin_r;
    Real cos_r;
    sincos( angle_r, &sin_r, &cos_r );  // do sincos here; latency. 

    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real hsigma( h * sigma );

    const Real num1( alpha * ( hsigma * sigma - hsigma * r * cos_r 
                               - ( r - sigma ) * cos_r ) 
                     + ( hsigma_p_1 + r * sigma * alphasq ) * sin_r );

    const Real num2( num_r0 );

    const Real den( r0 * alphasq * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( 2 * num1 * num2 / den );

    return result;
}


void 
FirstPassagePairGreensFunction::
createPsurvTable( RealVector& table, const Real r0 ) const
{
    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );

    table.clear();
    table.reserve( alphaTable_0.size() );

    std::transform( alphaTable_0.begin(), alphaTable_0.end(),
		    std::back_inserter( table ),
		    boost::bind( &FirstPassagePairGreensFunction::p_survival_i,
				 this, _1, r0 ) );
}


void 
FirstPassagePairGreensFunction::createNum_r0Table( RealVector& table,
						   const Real r0 ) const
{
    const RealVector& alphaTable_0( this->alphaTable[0] );

    table.clear();
    table.reserve( alphaTable_0.size() );

    std::transform( alphaTable_0.begin(), alphaTable_0.end(),
		    std::back_inserter( table ),
		    boost::bind( &FirstPassagePairGreensFunction::num_r0,
				 this, _1, r0 ) );
}

void 
FirstPassagePairGreensFunction::
createPleaveFactorTable( RealVector& table, const Real r0 ) const
{
    const RealVector& alphaTable_0( this->alphaTable[0] );

    table.clear();
    table.reserve( alphaTable_0.size() );

    std::transform( alphaTable_0.begin(), alphaTable_0.end(),
		    std::back_inserter( table ),
		    boost::bind( &FirstPassagePairGreensFunction::pleaveFactor,
				 this, _1, r0 ) );
}


void 
FirstPassagePairGreensFunction::
createPleavesTable( RealVector& table, const Real r0,
                    const RealVector& pleaveFactorTable ) const
{
    const RealVector& alphaTable_0( this->alphaTable[0] );

    assert( pleaveFactorTable.size() >= alphaTable_0.size() );

    table.clear();
    table.reserve( alphaTable_0.size() );

    for( unsigned int i( 0 ); i < alphaTable_0.size(); ++i )
    {
        const Real alpha( alphaTable_0[i] );
        table.push_back( p_leaves_i( alpha, r0, pleaveFactorTable[i] ) );
    }
}

void 
FirstPassagePairGreensFunction::
createPleaveaTable( RealVector& table, const Real r0,
                    const RealVector& pleaveFactorTable ) const
{
    const RealVector& alphaTable_0( this->alphaTable[0] );

    assert( pleaveFactorTable.size() >= alphaTable_0.size() );

    table.clear();
    table.reserve( alphaTable_0.size() );

    for( unsigned int i( 0 ); i < alphaTable_0.size(); ++i )
    {
        const Real alpha( alphaTable_0[i] );
        table.push_back( p_leavea_i( alpha, r0, pleaveFactorTable[i] ) );
    }
}


const Real 
FirstPassagePairGreensFunction::p_0_i_exp( const unsigned int i,
					   const Real t,
					   const Real r,
					   const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * p_0_i( alpha, r, r0 );
}


const Real 
FirstPassagePairGreensFunction::p_survival_i_exp( const unsigned int i,
						  const Real t,
						  const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * 
	p_survival_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::
p_survival_i_exp_table( const unsigned int i,
			const Real t,
			const Real r0,
			const RealVector& table ) const
{
    const Real alpha( this->getAlpha0( i ) );
//    printf("t %g\n",std::exp( - getD() * t * alpha * alpha ) * psurvTable[i] );
    return std::exp( - getD() * t * alpha * alpha ) * table[i];
}

const Real 
FirstPassagePairGreensFunction::
p_leave_i_exp_table( const unsigned int i,
                     const Real t,
                     const Real r0,
                     const RealVector& table ) const
{
    const Real alpha( this->getAlpha0( i ) );
//    printf("t %g\n",std::expm1( - getD() * t * alpha * alpha ) * psurvTable[i] );
    return expm1( - getD() * t * alpha * alpha ) * table[i];
}


const Real 
FirstPassagePairGreensFunction::dp_survival_i_exp( const unsigned int i,
						   const Real t,
						   const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * 
	dp_survival_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::leavea_i_exp( const unsigned int i,
					      const Real t,
					      const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * leavea_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::leaves_i_exp( const unsigned int i,
					      const Real t,
					      const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );

    return std::exp( - getD() * t * alpha * alpha ) * leaves_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::p_leavea_i_exp( const unsigned int i,
						const Real t,
						const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    const Real num_r0( this->num_r0( alpha, r0 ) ); 
    const Real den( this->p_survival_den( alpha, r0 ) ); 
    return exp( - getD() * t * alpha * alpha ) * 
        p_leavea_i( alpha, r0, num_r0 / den );
}

const Real 
FirstPassagePairGreensFunction::p_leaves_i_exp( const unsigned int i,
						const Real t,
						const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    const Real num_r0( this->num_r0( alpha, r0 ) ); 
    const Real den( this->p_survival_den( alpha, r0 ) ); 
    return exp( - getD() * t * alpha * alpha ) * 
        p_leaves_i( alpha, r0, num_r0 / den );
}

const Real 
FirstPassagePairGreensFunction::
p_int_r_i_exp( const unsigned int i,
	       const Real t,
	       const Real r,
	       const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );

    return std::exp( - getD() * t * alpha * alpha ) * 
	p_int_r_i( r, alpha, r0, num_r0( alpha, r0 ) );
}

const Real 
FirstPassagePairGreensFunction::
p_int_r_i_exp_table( const unsigned int i,
		     const Real t,
		     const Real r,
		     const Real r0,
		     const RealVector& num_r0Table ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * 
	p_int_r_i( r, alpha, r0, num_r0( alpha, r0 ) );//num_r0Table[i] );
}

const Real 
FirstPassagePairGreensFunction::p_0( const Real t,
				     const Real r,
				     const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_0_i_exp,
					this,
					_1, t, r, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}


const unsigned int
FirstPassagePairGreensFunction::guess_maxi( const Real t ) const
{
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real a( geta() );

    const Real Dt( D * t );
    const Real max_alpha( sqrt( - log( this->TOLERANCE ) / Dt ) );
    
    return static_cast<unsigned int>( max_alpha * ( a - sigma ) / M_PI ) + 1;
}


const Real 
FirstPassagePairGreensFunction::p_survival( const Real t,
					    const Real r0 ) const
{
/*
    unsigned int maxi( guess_maxi( t ) );
    printf("maxi %d\n",maxi );
    if( maxi > this->MAX_ALPHA_SEQ )
    {

    }
*/

    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_survival_i_exp, 
					this,
					_1, t, r0 ), 
			   this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::
p_survival_table( const Real t,
		  const Real r0,
		  const RealVector& psurvTable ) const
{
    const Real p( funcSum_all( boost::bind( &FirstPassagePairGreensFunction::
                                            p_survival_i_exp_table, 
                                            this,
                                            _1, t, r0, psurvTable ),
                               psurvTable.size() ) );

    return p;
}

const Real 
FirstPassagePairGreensFunction::
p_leave_table( const Real t,
               const Real r0,
               const RealVector& table ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_leave_i_exp_table, 
					this,
					_1, t, r0, table ),
			   table.size() ) );

    return p;
}


const Real 
FirstPassagePairGreensFunction::dp_survival( const Real t,
					     const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					dp_survival_i_exp, 
					this,
					_1, t, r0 ),
                                        this->MAX_ALPHA_SEQ ) );
    return p;
}


const Real 
FirstPassagePairGreensFunction::leaves( const Real t,
					const Real r0 ) const
{
/*
    this->updateAlphaTable0( t );
    const Real p( funcSum_all_accel( boost::bind( &FirstPassagePairGreensFunction::
                                                  leaves_i_exp,
                                                  this,
                                                  _1, t, r0 ),
                                     this->getAlphaTable( 0 ).size() ) );*/
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
                                            leaves_i_exp,
                                            this,
                                            _1, t, r0 ),
                               this->MAX_ALPHA_SEQ ) );

    return p;
}

const Real 
FirstPassagePairGreensFunction::leavea( const Real t,
					const Real r0 ) const
{
/*
    this->updateAlphaTable0( t );
    const Real p( funcSum_all_accel( boost::bind( &FirstPassagePairGreensFunction::
                                            leavea_i_exp,
                                            this,
                                            _1, t, r0 ),
                                     this->getAlphaTable( 0 ).size() ) );
*/

    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
                                            leavea_i_exp,
                                            this,
                                            _1, t, r0 ),
                               this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::p_leaves( const Real t,
					  const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_leaves_i_exp,
					this,
					_1, t, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}


const Real 
FirstPassagePairGreensFunction::p_leavea( const Real t,
					  const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_leavea_i_exp,
					this,
					_1, t, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::p_int_r( const Real r,
					 const Real t,
					 const Real r0 ) const

{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_int_r_i_exp,
					this,
					_1, t, r, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::
p_int_r_table( const Real r,
	       const Real t,
	       const Real r0,
	       const RealVector& num_r0Table ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_int_r_i_exp_table,
					this,
					_1, t, r, r0, num_r0Table ),
			   num_r0Table.size() ) );
    return p;
}



const Real
FirstPassagePairGreensFunction::p_survival_F( const Real t,
					      const p_survival_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const RealVector& table( params->table );
    const Real rnd( params->rnd );

    return rnd - gf->p_survival_table( t, r0, table );
}

const Real
FirstPassagePairGreensFunction::p_leave_F( const Real t,
                                            const p_survival_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const RealVector& table( params->table );
    const Real rnd( params->rnd );
//    printf("p_leaves_F %g %g %g %g\n", t, rnd, gf->p_survival_table( t, r0, psurvTable ),
    //         gf->p_survival_table( t, r0, psurvTable ) - rnd );
    return - gf->p_leave_table( t, r0, table ) - rnd;
}


const Real
FirstPassagePairGreensFunction::p_int_r_F( const Real r,
					   const p_int_r_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real r0( params->r0 );
//    const RealVector& num_r0Table( params->num_r0Table );
    const Real rnd( params->rnd );

//    return gf->p_int_r_table( r, t, r0, num_r0Table ) - rnd;
    return gf->p_int_r( r, t, r0 ) - rnd;
}


const Real FirstPassagePairGreensFunction::drawTime( const Real rnd, 
						     const Real r0 ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 <= a );

    if( r0 == a || a == sigma )
    {
	return 0.0;
    }

    Real t_guess;

    //const Real dist( std::min( a - r0, r0 - sigma ) );
    const Real dist( a - r0 );  // this is a good guess, actually.
    t_guess = dist * dist / ( 6.0 * D );
    t_guess *= .1;

    const Real minT( std::min( sigma * sigma / D * this->MIN_T_FACTOR,
                               t_guess * 1e-5 ) );

    RealVector psurvTable;

    p_survival_params params = { this, r0, psurvTable, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_F ),
	    &params 
	};

    this->updateAlphaTable0( t_guess );
    this->createPsurvTable( psurvTable, r0 );

    Real low( t_guess );
    Real high( t_guess );

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    const Real value( GSL_FN_EVAL( &F, t_guess ) );
    if( value < 0.0 )
    {
        high *= 10;
        while( 1 )
        {
            const Real high_value( GSL_FN_EVAL( &F, high ) );
            
            if( high_value >= 0.0 )
            {
                break;
            }

            if( fabs( high ) >= 1e10 )
            {
                std::cerr << "Couldn't adjust high. F(" << high <<
                    ") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
                    ", " << dump() << std::endl;
                throw std::exception();
            }

//            printf( "drawTime: adjusting high: %g F = %g\n", high, high_value );
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev( value );
        low *= .1;

        while( 1 )
        {
            this->updateAlphaTable0( low );
            this->createPsurvTable( psurvTable, r0 );
            
            const Real low_value( GSL_FN_EVAL( &F, low ) );
            
            if( low_value <= 0.0 )
            {
                break;
            }
            
            // FIXME: 
            if( fabs( low ) <= minT || 
                fabs( low_value - low_value_prev ) < TOLERANCE ) 
            {
                std::cerr << "Couldn't adjust low.  Returning minT (= "
                          << minT << "); F(" << low <<
                    ") = " << GSL_FN_EVAL( &F, low ) << "; r0 = " << r0 << ", "
                          << dump() << std::endl;
                return minT;
            }
            low_value_prev = low_value;

//            printf( "drawTime: adjusting low: %g, F = %g\n", low, low_value );
            low *= .1;
        }
    }

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real t( findRoot( F, solver, low, high, this->TOLERANCE, 
                            "drawTime" ) );

    gsl_root_fsolver_free( solver );

    return t;
}

const EventType
FirstPassagePairGreensFunction::drawEventType( const Real rnd, 
					       const Real r0,
					       const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real kf( this->getkf() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    if( kf == 0 )
    {
        return ESCAPE;
    }

    const Real reaction( leaves( t, r0 ) * 4.0 * M_PI * sigma * sigma );
    const Real escape( leavea( t, r0 ) * 4.0 * M_PI * a * a );
    //const Real den( dp_survival( t, r0 ) );

    //const Real value( reaction / den );
    const Real value( reaction / ( reaction + escape ) );
    //printf("et %g %g %g %g\n", value, reaction, escape, den );

    //assert( value >= - 1e-7 && value <= 1.0 + 1e-7 );

    if( rnd <= value )  
    {
	return REACTION;   // leaves
    }
    else 
    {
	return ESCAPE;     // leavea
    }
}


const boost::tuple<Real,EventType>
FirstPassagePairGreensFunction::drawTime2( const Real rnd1, 
                                           const Real rnd2,
                                           const Real r0 ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );
    const Real kf( this->getkf() );

    THROW_UNLESS( std::invalid_argument, rnd1 < 1.0 && rnd1 >= 0.0 );
    THROW_UNLESS( std::invalid_argument, rnd2 < 1.0 && rnd2 >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 <= a );

    if( r0 == a || a == sigma )
    {
	return boost::make_tuple( 0.0, ESCAPE );
    }

    RealVector pleaveFactorTable;
    RealVector pleavesTable;
    RealVector pleaveaTable;

    const Real pleavea_inf( p_leavea( INFINITY, r0 ) );
    p_survival_params params_a = { this, r0, pleaveaTable, rnd1 * pleavea_inf };
    gsl_function Fa = 
	{
	    reinterpret_cast<typeof(Fa.function)>( &p_leave_F ),
	    &params_a 
	};

    Real t_guess;

    if( kf != 0.0 )
    {
        const Real dist( std::min( a - r0, r0 - sigma ) );
        t_guess = dist * dist / ( 6.0 * D );
    }
    else // kf == 0.0
    {
        const Real dist( a - r0 );
        t_guess = dist * dist / ( 6.0 * D );
    }

    const Real minT( std::min( sigma * sigma / D * this->MIN_T_FACTOR,
                               t_guess * 1e-2 ) );

    /*
    if( kf == 0.0 )
    {
        gsl_root_fsolver_free( solver );
        return boost::make_tuple( findRoot( Fa, solver, low, high, 
                                            this->TOLERANCE,
                                            "drawTime2: escape" ),
                                  ESCAPE );
    }
    */


    const Real pleaves_inf( p_leaves( INFINITY, r0 ) );
    p_survival_params params_s = { this, r0, pleavesTable, rnd2 * pleaves_inf };
    gsl_function Fs = 
	{
	    reinterpret_cast<typeof(Fs.function)>( &p_leave_F ),
	    &params_s
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    this->updateAlphaTable0( t_guess * 1e-3 );
    this->createPleaveFactorTable( pleaveFactorTable, r0 );
    this->createPleavesTable( pleavesTable, r0, pleaveFactorTable );
    this->createPleaveaTable( pleaveaTable, r0, pleaveFactorTable );

    Real value_s( GSL_FN_EVAL( &Fs, t_guess ) );
    Real value_a( GSL_FN_EVAL( &Fa, t_guess ) );
    Real high( t_guess );
    Real low( t_guess );

    if( value_s <= 0.0 )
    {
        if( value_a <= 0.0 )
        {
            // need to adjust high for both
            while( 1 )
            {
                high *= 10;
                if( fabs( high ) >= 1e10 )
                {
                    std::cerr << "Couldn't adjust high (reaction). F(" 
                              << high << ") = " << GSL_FN_EVAL( &Fs, high ) 
                              << "; r0 = " << r0 
                              << ", " << dump() << std::endl;
                    gsl_root_fsolver_free( solver );
                    throw std::exception();
                }

                value_s = GSL_FN_EVAL( &Fs, high );
                value_a = GSL_FN_EVAL( &Fa, high );

                printf( "drawTime2: adjusting high: %g, Fs = %g, Fa= %g\n", 
                        high, value_s, value_a );

                if( value_s >= 0.0 && value_a >= 0.0 )
                {
                    // both straddle
                    break;
                }

                if( value_s * value_a <= 0.0 )
                {
                    if( value_s >= 0.0 )
                    {
                        boost::tuple<Real,EventType> 
                            ret(  boost::make_tuple( findRoot( Fs, solver, 
                                                               low, high, 
                                                               this->TOLERANCE,
                                                               "drawTime2: s" ),
                                                     REACTION ) );
                        gsl_root_fsolver_free( solver );
                        return ret;
                    }
                    else // value_a >= 0.0
                    {
                        boost::tuple<Real,EventType>
                            ret( boost::make_tuple( findRoot( Fa, solver, 
                                                              low, high, 
                                                               this->TOLERANCE,
                                                              "drawTime2: a" ),
                                                    ESCAPE ) );
                        gsl_root_fsolver_free( solver );
                        return ret;
                    }
                }
          }
        }
        else  // value_s <= 0.0, value_a > 0.0
        {
            // adjust low for only a

            while( 1 )
            {
                low *= .1;

                if( fabs( low ) <= minT )//|| 
                    //fabs( low_value - low_value_prev ) < TOLERANCE ) 
                {
                    std::cerr << "Couldn't adjust low.  Returning minT (= "
                              << minT << "); F(" << low 
                              << ") = " << GSL_FN_EVAL( &Fa, low ) 
                              << "; r0 = " << r0 << ", "
                              << dump() << std::endl;

                    gsl_root_fsolver_free( solver );
                    return boost::make_tuple( minT, ESCAPE ); // FIXME
                }

                this->updateAlphaTable0( low );
                this->createPleaveFactorTable( pleaveFactorTable, r0 );
                this->createPleaveaTable( pleaveaTable, r0, pleaveFactorTable );

                const Real value_a( GSL_FN_EVAL( &Fa, low ) );

                printf( "drawTime2: adjusting low: %g, Fa = %g\n", 
                        low, value_a );
                printf("%g %g\n",value_a, value_s );

                if( value_a <= 0.0 )
                {
                    boost::tuple<Real,EventType> 
                        ret( boost::make_tuple( findRoot( Fa, solver, 
                                                          low, high, 
                                                          this->TOLERANCE,
                                                          "drawTime2: s" ),
                                                ESCAPE ) );
                    gsl_root_fsolver_free( solver );
                    return ret;
                }
            }
        }
    }
    else // value_s > 0.0
    {
        if( value_a > 0.0 )
        {
          // adjust low for both
            while( 1 )
            {
                low *= .1;

                if( fabs( low ) <= minT )//|| 
                    //fabs( low_value - low_value_prev ) < TOLERANCE ) 
                {
                    std::cerr << "Couldn't adjust low.  Returning minT (= "
                              << minT << "); F(" << low 
                              << ") = " << value_s
                              << "; r0 = " << r0 << ", "
                              << dump() << std::endl;
                    gsl_root_fsolver_free( solver );
                    return boost::make_tuple( minT, 
                                              value_s > value_a ?
                                              REACTION : ESCAPE ); // FIXME
                }

                this->updateAlphaTable0( low );
                this->createPleaveFactorTable( pleaveFactorTable, r0 );
                this->createPleavesTable( pleavesTable, r0, pleaveFactorTable );
                this->createPleaveaTable( pleaveaTable, r0, pleaveFactorTable );

                value_s = GSL_FN_EVAL( &Fs, low );
                value_a = GSL_FN_EVAL( &Fa, low );

                printf( "drawTime2: adjusting low: %g, Fs = %g, Fa = %g\n", 
                        low, value_s, value_a );

                if( value_s <= 0.0 && value_a <= 0.0 )
                {
                    break;
                }

                if( value_s * value_a <= 0.0 )
                {
                    if( value_s <= 0.0 )
                    {
                        boost::tuple<Real,EventType> 
                            ret( boost::make_tuple( findRoot( Fs, solver, 
                                                              low, high, 
                                                               this->TOLERANCE,
                                                              "drawTime2: s" ),
                                                    REACTION ) );
                        gsl_root_fsolver_free( solver );
                        return ret;
                    }
                    else // value_a <= 0.0
                    {
                        boost::tuple<Real,EventType> 
                            ret( boost::make_tuple( findRoot( Fa, solver, 
                                                              low, high, 
                                                               this->TOLERANCE,
                                                              "drawTime2: a" ),
                                                    ESCAPE ) );
                        gsl_root_fsolver_free( solver );
                        return ret;
                    }
                }
            }
        }
        else // value_s > 0.0, value_a <= 0.0
        {
            // adjust low for only s

            while( 1 )
            {
                low *= .1;

                this->updateAlphaTable0( low );
                this->createPleaveFactorTable( pleaveFactorTable, r0 );
                this->createPleavesTable( pleavesTable, r0, pleaveFactorTable );

                const Real value_s( GSL_FN_EVAL( &Fs, low ) );

                printf( "drawTime2: adjusting low: %g, Fs = %g\n", 
                        low, value_s );
                    printf("%g %g\n",value_a, value_s );

                if( value_s <= 0.0 )
                {
                    boost::tuple<Real,EventType> 
                        ret( boost::make_tuple( findRoot( Fs, solver, 
                                                          low, high, 
                                                          this->TOLERANCE,
                                                          "drawTime2: s" ),
                                                REACTION ) );
                    gsl_root_fsolver_free( solver );
                    return ret;
                }

                if( fabs( low ) <= minT )//|| 
                    //fabs( low_value - low_value_prev ) < TOLERANCE ) 
                {
                    std::cerr << "Couldn't adjust low.  Returning minT (= "
                              << minT << "); F(" << low 
                              << ") = " << GSL_FN_EVAL( &Fs, low ) 
                              << "; r0 = " << r0 << ", "
                              << dump() << std::endl;

                    gsl_root_fsolver_free( solver );
                    return boost::make_tuple( minT, REACTION );
                }
            }
        }
    }

    const Real t_escape( findRoot( Fa, solver, low, high, this->TOLERANCE,
                                   "drawTime2: escape" ) );
    const Real t_reaction( findRoot( Fs, solver, low, high, this->TOLERANCE,
                                     "drawTime2: reaction" ) );

    printf("e %g r %g\n", t_escape, t_reaction);

    boost::tuple<Real,EventType> ret;
    if( t_escape <= t_reaction )
    {
        ret = boost::make_tuple( t_escape, ESCAPE );
    }
    else
    {
        ret = boost::make_tuple( t_reaction, REACTION );
    }

    gsl_root_fsolver_free( solver );

    return ret;
}


#if 0
const boost::tuple<Real,EventType>
FirstPassagePairGreensFunction::drawTime2( const Real rnd1, 
                                           const Real rnd2,
                                           const Real r0 ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );
    const Real kf( this->getkf() );

    THROW_UNLESS( std::invalid_argument, rnd1 < 1.0 && rnd1 >= 0.0 );
    THROW_UNLESS( std::invalid_argument, rnd2 < 1.0 && rnd2 >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 <= a );

    if( r0 == a || a == sigma )
    {
	return boost::make_tuple( 0.0, ESCAPE );
    }

    //const Real minT( std::min( sigma * sigma / D * this->MIN_T_FACTOR,
//                               t_guess * 1e-2 ) );

    RealVector pleaveFactorTable;
    RealVector pleavesTable;
    RealVector pleaveaTable;


    // INCOMPLETE: normalize
    p_survival_params params_a = { this, r0, pleaveaTable, rnd1 };
    gsl_function Fa = 
	{
	    reinterpret_cast<typeof(Fa.function)>( &p_leave_F ),
	    &params_a 
	};

    // INCOMPLETE: normalize
    p_survival_params params_s = { this, r0, pleavesTable, rnd2 };
    gsl_function Fs = 
	{
	    reinterpret_cast<typeof(Fs.function)>( &p_leave_F ),
	    &params_s
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real dist_a( a - r0 );
    const Real t_guess_a( dist_a * dist_a / ( 6.0 * D ) );
    Real t_escape;

    if( kf == 0.0 )
    {
        this->updateAlphaTable0( t_guess_a );
        this->createPleaveFactorTable( pleaveFactorTable, r0 );
        this->createPleaveaTable( pleaveaTable, r0, pleaveFactorTable );

        t_escape = drawPleavea( Fa, solver, r0, t_guess_a,
                                pleaveFactorTable,
                                pleaveaTable );

        gsl_root_fsolver_free( solver );
        return boost::make_tuple( t_escape, ESCAPE );
    }

    const Real dist_s( r0 - sigma );
    const Real t_guess_s( dist_s * dist_s / ( 6.0 * D ) );
    Real t_reaction;

    if( dist_a <= dist_s )
    {
        this->updateAlphaTable0( t_guess_a );
        this->createPleaveFactorTable( pleaveFactorTable, r0 );
        this->createPleaveaTable( pleaveaTable, r0, pleaveFactorTable );

        t_escape = drawPleavea( Fa, solver, r0, t_guess_a, 
                                pleaveFactorTable,
                                pleaveaTable );

        this->updateAlphaTable0( t_guess_s );
        this->createPleaveFactorTable( pleaveFactorTable, r0 );
        this->createPleavesTable( pleavesTable, r0, pleaveFactorTable );

        if( t_escape < t_guess_s )
        {
            const Real value_s( GSL_FN_EVAL( &Fs, t_guess_s ) );
            if( value_s <= 0.0 )
            {
                gsl_root_fsolver_free( solver );
                return boost::make_tuple( t_escape, ESCAPE );
            }
        }

        t_reaction = drawPleaves( Fs, solver, r0, t_guess_s,
                                  pleaveFactorTable,
                                  pleavesTable );
    }
    else
    {
        this->updateAlphaTable0( t_guess_s );
        this->createPleaveFactorTable( pleaveFactorTable, r0 );
        this->createPleavesTable( pleavesTable, r0, pleaveFactorTable );

        t_reaction = drawPleaves( Fs, solver, r0, t_guess_s,
                                  pleaveFactorTable,
                                  pleavesTable );

        this->updateAlphaTable0( t_guess_a );
        this->createPleaveFactorTable( pleaveFactorTable, r0 );
        this->createPleaveaTable( pleaveaTable, r0, pleaveFactorTable );

        if( t_reaction < t_guess_a )
        {
            const Real value_a( GSL_FN_EVAL( &Fa, t_guess_a ) );
            if( value_a <= 0.0 )
            {
                gsl_root_fsolver_free( solver );
                return boost::make_tuple( t_reaction, REACTION );
            }
        }

        t_escape = drawPleavea( Fa, solver, r0, t_guess_a,
                                pleaveFactorTable,
                                pleaveaTable );
    }

    printf("e %g r %g\n", t_escape, t_reaction);

    boost::tuple<Real,EventType> ret;
    if( t_escape <= t_reaction )
    {
        ret = boost::make_tuple( t_escape, ESCAPE );
    }
    else
    {
        ret = boost::make_tuple( t_reaction, REACTION );
    }

    gsl_root_fsolver_free( solver );

    return ret;
}
#endif

const Real 
FirstPassagePairGreensFunction::drawPleavea( gsl_function& F,
                                             gsl_root_fsolver* solver,
                                             const Real r0,
                                             const Real t_guess,
                                             RealVector& pleaveFactorTable,
                                             RealVector& pleaveaTable ) const
{
    const Real minT( 1e-12 );

    Real low( t_guess );
    Real high( t_guess );

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    const Real value( GSL_FN_EVAL( &F, t_guess ) );
    if( value < 0.0 )
    {
        high *= 10;
        while( 1 )
        {
            const Real high_value( GSL_FN_EVAL( &F, high ) );
            
            if( high_value >= 0.0 )
            {
                break;
            }

            if( fabs( high ) >= 1e10 )
            {
                std::cerr << "Couldn't adjust high. Fa(" << high <<
                    ") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
                    ", " << dump() << std::endl;
                throw std::exception();
            }

            printf( "drawTime2: adjusting high: %g Fa = %g\n", 
                    high, high_value );
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev( value );
        low *= .1;

        while( 1 )
        {
            this->updateAlphaTable0( low );
            this->createPleaveFactorTable( pleaveFactorTable, r0 );
            this->createPleaveaTable( pleaveaTable, r0, pleaveFactorTable );

            
            const Real low_value( GSL_FN_EVAL( &F, low ) );
            
            if( low_value <= 0.0 )
            {
                break;
            }
            
            // FIXME: 
            if( fabs( low ) <= minT || 
                fabs( low_value - low_value_prev ) < TOLERANCE ) 
            {
                std::cerr << "Couldn't adjust low.  Returning minT (= "
                          << minT << "); Fa(" << low <<
                    ") = " << GSL_FN_EVAL( &F, low ) << "; r0 = " << r0 << ", "
                          << dump() << std::endl;
                return minT;
            }
            low_value_prev = low_value;

            printf( "drawTime2: adjusting low: %g, Fa = %g\n", low, low_value );
            low *= .1;
        }
    }

    const Real t( findRoot( F, solver, low, high,  
                            this->TOLERANCE, "drawTime2: a" ) );

    return t;
}


const Real 
FirstPassagePairGreensFunction::drawPleaves( gsl_function& F,
                                             gsl_root_fsolver* solver,
                                             const Real r0,
                                             const Real t_guess,
                                             RealVector& pleaveFactorTable,
                                             RealVector& pleavesTable ) const
{
    const Real minT( 1e-12 );

    Real low( t_guess );
    Real high( t_guess );

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    const Real value( GSL_FN_EVAL( &F, t_guess ) );
    if( value < 0.0 )
    {
        high *= 10;
        while( 1 )
        {
            const Real high_value( GSL_FN_EVAL( &F, high ) );
            
            if( high_value >= 0.0 )
            {
                break;
            }

            if( fabs( high ) >= 1e10 )
            {
                std::cerr << "Couldn't adjust high. Fs(" << high <<
                    ") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
                    ", " << dump() << std::endl;
                throw std::exception();
            }

            printf( "drawTime2: adjusting high: %g Fs = %g\n", 
                    high, high_value );
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev( value );
        low *= .1;

        while( 1 )
        {
            this->updateAlphaTable0( low );
            this->createPleaveFactorTable( pleaveFactorTable, r0 );
            this->createPleavesTable( pleavesTable, r0, pleaveFactorTable );
            
            const Real low_value( GSL_FN_EVAL( &F, low ) );
            
            if( low_value <= 0.0 )
            {
                break;
            }
            
            // FIXME: 
            if( fabs( low ) <= minT )//|| 
//                fabs( low_value - low_value_prev ) < TOLERANCE ) 
            {
                std::cerr << "Couldn't adjust low.  Returning minT (= "
                          << minT << "); Fs(" << low <<
                    ") = " << GSL_FN_EVAL( &F, low ) << "; r0 = " << r0 << ", "
                          << dump() << std::endl;
                return minT;
            }
            low_value_prev = low_value;

            printf( "drawTime2: adjusting low: %g, Fs = %g\n", low, low_value );
            low *= .1;
        }
    }

    const Real t( findRoot( F, solver, low, high, this->TOLERANCE,
                            "drawTime2: s" ) );

    return t;
}




const Real FirstPassagePairGreensFunction::drawR( const Real rnd, 
						  const Real r0, 
						  const Real t ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );

    if( t == 0.0 )
    {
	return r0;
    }

    const Real psurv( p_survival( t, r0 ) );

//    RealVector num_r0Table;
//    createNum_r0Table( num_r0Table, r0 );

    p_int_r_params params = { this, t, r0, /*num_r0Table,*/ rnd * psurv };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};


    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.

    Real low( r0 );
    Real high( r0 );

    const Real sqrt6Dt( sqrt( 6.0 * D * t ) );
    if( GSL_FN_EVAL( &F, r0 ) < 0.0 )
    {
        // low = r0
        unsigned int H( 3 );

        while( true )
        {
            high = r0 + H * sqrt6Dt;
            if( high > a )
            {
                if( GSL_FN_EVAL( &F, a ) < 0.0 )
                {
                    printf( "drawR: p_int_r( a ) < 0.0. returning a.\n" );
                    return a;
                }

                high = a;
                break;
            }

            const Real value( GSL_FN_EVAL( &F, high ) );
            if( value > 0.0 )
            {
                break;
            }

            ++H;
        }

    }
    else
    {
        // high = r0
        unsigned int H( 3 );

        while( true )
        {
            low = r0 - H * sqrt6Dt;
            if( low < sigma )
            {
                if( GSL_FN_EVAL( &F, sigma ) > 0.0 )
                {
                    printf( "drawR: p_int_r( sigma ) > 0.0. "
                            "returning sigma.\n" );
                    return sigma;
                }

                low = sigma;
                break;
            }

            const Real value( GSL_FN_EVAL( &F, low ) );
            if( value < 0.0 )
            {
                break;
            }

            ++H;
        }
    }


    // root finding by iteration.

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	low = gsl_root_fsolver_x_lower( solver );
	high = gsl_root_fsolver_x_upper( solver );
	const int status( gsl_root_test_interval( low, high, 1e-15,
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawR: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    //printf("%d\n", i );


    const Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return r;
}




const Real FirstPassagePairGreensFunction::p_n_alpha( const unsigned int i,
						      const unsigned int n,
						      const Real r,
						      const Real r0, 
						      const Real t ) const
//                                                      const RealVector& jas1Table ) const
{
    const Real sigma( this->getSigma() );
    const Real h( this->geth() );
    const Real a( this->geta() );

    const Real mDt( - this->getD() * t );

    const Real alpha( this->getAlpha( n, i ) );
    const Real alphasq( alpha * alpha );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );
    const Real hSigma_m_n( hSigma - realn );

    const Real term1( alphasq * alphasq * exp( mDt * alphasq ) );

    // GSL
    const Real jas1( gsl_sf_bessel_jl( n,   sigmaAlpha ) );
    const Real yas1( gsl_sf_bessel_yl( n,   sigmaAlpha ) );
    const Real jas2( gsl_sf_bessel_jl( n+1, sigmaAlpha ) );
    const Real yas2( gsl_sf_bessel_yl( n+1, sigmaAlpha ) );
    const Real jaa(  gsl_sf_bessel_jl( n,   aAlpha ) );
    const Real yaa(  gsl_sf_bessel_yl( n,   aAlpha ) );
    const Real jar(  gsl_sf_bessel_jl( n,   r * alpha ) );
    const Real yar(  gsl_sf_bessel_yl( n,   r * alpha ) );
    const Real jar0( gsl_sf_bessel_jl( n,   r0 * alpha ) );
    const Real yar0( gsl_sf_bessel_yl( n,   r0 * alpha ) );

    const Real J( hSigma_m_n * jas1 + sigmaAlpha * jas2 );
    const Real Y( hSigma_m_n * yas1 + sigmaAlpha * yas2 );
    const Real falpha_r( - J * yar + Y * jar );
    const Real falpha_r0( - J * yar0 + Y * jar0 );

    const Real num( falpha_r * falpha_r0 );

    const Real E1( a * ( realn + realn * realn - 
                         sigma * ( h + h * h * sigma + sigma * alphasq ) ) );

    const Real E2( sigma * ( J * J + Y * Y ) / ( jaa * jaa + yaa * yaa ) );
    //const Real E2( sigma * J * J / ( jaa * jaa ) );

    const Real den( E1 + E2 );

    const Real result( term1 * num / den );

    return result;
}



const Real 
FirstPassagePairGreensFunction::p_n( const Integer n,
				     const Real r,
				     const Real r0, 
				     const Real t ) const
{


    /*
    const Real sigma( getSigma() );

    RealVector& alphaTable( getAlphaTable( n ) );
    
    Real realn( 0.5 + n );
    RealVector jas1Table( alphaTable );

    for( RealVector::iterator i( alphaTable.begin() );  i !=alphaTable.end();
         ++i )
    {
        *i = *i * sigma;
    }

    gsl_sf_bessel_sequence_Jnu_e( realn, GSL_PREC_SINGLE, 
                                  jas1Table.size(), &jas1Table[0] );
    */

    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_n_alpha,
					this,
					_1, n, r, r0, t ),// jas1Table ),
			   this->MAX_ALPHA_SEQ ) );


    return p;
}

void
FirstPassagePairGreensFunction::makep_nTable( RealVector& p_nTable,
					      const Real r, 
					      const Real r0, 
					      const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    p_nTable.clear();

    const Real factor( a * sigma / ( M_PI * 2 ) );

    const Real p_0( this->p_n( 0, r, r0, t ) * factor );
    p_nTable.push_back( p_0 );

    const Real threshold( fabs( this->TOLERANCE * p_0 * 1e-2 ) );

    Real p_n_prev_abs( fabs( p_0 ) );
    unsigned int n( 1 );
    while( true )
    {
	Real p_n( this->p_n( n, r, r0, t ) * factor );

	if( ! std::isfinite( p_n ) )
	{
	    std::cerr << "makep_nTable: invalid value; " <<
		p_n << "( n= " << n << ")." << std::endl;
//	    p_n = 0.0;
	    break;
	}
	//printf("%d p_n %g\n", n, p_n );

	p_nTable.push_back( p_n );
        
        // std::cerr << n << " " << p_n << " " << threshold << std::endl;

	const Real p_n_abs( fabs( p_n ) );
	// truncate when converged enough.
	if( p_n_abs < threshold &&
            p_n_prev_abs < threshold &&
	    p_n_abs <= p_n_prev_abs &&
            n >= 4 )  // corresponds to funcSum()'s CONVERGENCE_CHECK
	{
	    break;
        }
	
	++n;

	if( n >= this->MAX_ORDER )
	{
	    // std::cerr << "p_n didn't converge." << std::endl;
	    break;
	}
	
	p_n_prev_abs = p_n_abs;
    }

}


const Real 
FirstPassagePairGreensFunction::dp_n_alpha_at_a( const unsigned int i,
						 const unsigned int n,
						 const Real r0, 
						 const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real h( this->geth() );
    const Real a( this->geta() );

    const Real mDt( - this->getD() * t );
    const Real alpha( this->getAlpha( n, i ) );

    const Real alphasq( alpha * alpha );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );
    const Real hSigma_m_n( hSigma - realn );

    const Real term1( alphasq * alpha * exp( mDt * alphasq ) );

    // GSL
    const Real jas1( gsl_sf_bessel_jl( n,   sigmaAlpha ) );
    const Real yas1( gsl_sf_bessel_yl( n,   sigmaAlpha ) );
    const Real jas2( gsl_sf_bessel_jl( n+1, sigmaAlpha ) );
    const Real yas2( gsl_sf_bessel_yl( n+1, sigmaAlpha ) );
    const Real jaa1( gsl_sf_bessel_jl( n,   aAlpha ) );
    const Real yaa1( gsl_sf_bessel_yl( n,   aAlpha ) );
    const Real jar0( gsl_sf_bessel_jl( n,   r0 * alpha ) );
    const Real yar0( gsl_sf_bessel_yl( n,   r0 * alpha ) );

    const Real J( hSigma_m_n * jas1 + sigmaAlpha * jas2 );
    const Real Y( hSigma_m_n * yas1 + sigmaAlpha * yas2 );

    const Real falpha_r0( - J * yar0 + Y * jar0 );

    /*
    const Real num1( - J * jaa1 );
    const Real den1( a * ( realn + realn * realn - 
                           sigma * ( h + h * h * sigma + sigma * alphasq ) )
                     * jaa1 * jaa1 );
    const Real den2( sigma * J * J );
    const Real den( den1 + den2 );
    */

// Alternative form;

     const Real num1( - J );
     const Real den1( a * ( realn + realn * realn - 
                            sigma * ( h + h * h * sigma + sigma * alphasq ) ) );
     const Real den2( sigma * ( J * J + Y * Y ) / 
 		   ( jaa1 * jaa1 + yaa1 * yaa1 ) );
     const Real den( ( den1 + den2 ) * jaa1 );


     const Real num( num1 * falpha_r0 );


//    printf("f n1 n2 d %g %g %g %g\n",falpha_r0, num1, num2, den );

    const Real result( term1 * num / den );

    return result;
}

const Real 
FirstPassagePairGreensFunction::dp_n_at_a( const Integer n,
					   const Real r0, 
					   const Real t ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					dp_n_alpha_at_a,
					this,
					_1, n, r0, t ),
			   this->MAX_ALPHA_SEQ ) );

    return p;
}


void
FirstPassagePairGreensFunction::makedp_n_at_aTable( RealVector& p_nTable,
						    const Real r0, 
						    const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    p_nTable.clear();

    const Real factor( getD() * sigma / ( 2.0 * a * M_PI ) );

    const Real p_0( this->dp_n_at_a( 0, r0, t ) * factor );
    p_nTable.push_back( p_0 );

    const Real threshold( fabs( this->TOLERANCE * p_0 * 1e-2 ) );
    //printf("p_0 %g\n",p_0 );
    Real p_n_prev_abs( fabs( p_0 ) );
    unsigned int n( 1 );
    while( true )
    {
	Real p_n( this->dp_n_at_a( n, r0, t ) * factor );

	if( ! std::isfinite( p_n ) )
	{
	    std::cerr << "makedp_n_at_aTable: invalid value; " <<
		p_n << "( n= " << n << ")." << std::endl;
//	    p_n = 0.0;
	    break;
	}
	//printf("p_n %g\n",p_n );

	p_nTable.push_back( p_n );

	const Real p_n_abs( fabs( p_n ) );
	// truncate when converged enough.
	if( p_n_abs < threshold &&
            p_n_prev_abs < threshold &&
	    p_n_abs < p_n_prev_abs )
	{
	    break;
	}
	
	++n;

	if( n >= this->MAX_ORDER )
	{
	    //std::cerr << "dp_n_at_a didn't converge." << std::endl;
	    break;
	}
	
	p_n_prev_abs = p_n_abs;
    }

}

const Real 
FirstPassagePairGreensFunction::p_theta( const Real theta,
					 const Real r, 
					 const Real r0, 
					 const Real t ) const 
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in ( sigma, a );  not defined at r == sigma and r == a.
	THROW_UNLESS( std::invalid_argument, r >= sigma && r < a );
	THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makep_nTable( p_nTable, r, r0, t );

    const Real p( p_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}


const Real 
FirstPassagePairGreensFunction::dp_theta( const Real theta,
                                          const Real r, 
                                          const Real r0, 
                                          const Real t ) const 
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );

        // r \in [ sigma, a ]  ;  unlike p_theta,
        // defined at r == sigma and r == a.
	THROW_UNLESS( std::invalid_argument, r >= sigma && r <= a );
	THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makedp_n_at_aTable( p_nTable, r0, t );

    const Real p( p_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}


const Real 
FirstPassagePairGreensFunction::
p_theta_n( const unsigned int n,
	   const RealVector& p_nTable, const RealVector& lgndTable ) const
{
    return p_nTable[n] * lgndTable[n] * ( 2 * n + 1 );
}

const Real
FirstPassagePairGreensFunction::
p_theta_table( const Real theta,
	       const Real r, 
	       const Real r0, 
	       const Real t, 
	       const RealVector& p_nTable ) const
{
    const unsigned int tableSize( p_nTable.size() );

    Real sin_theta;
    Real cos_theta;
    sincos( theta, &sin_theta, &cos_theta );

    RealVector lgndTable( tableSize );
    gsl_sf_legendre_Pl_array( tableSize-1, cos_theta, &lgndTable[0] );

    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_theta_n,
					this,
					_1, p_nTable, lgndTable ),
			   tableSize ) );

    return p * sin_theta;
}


void
FirstPassagePairGreensFunction::
make_p_thetaTable( RealVector& pTable,
		   const Real r, 
		   const Real r0, 
		   const Real t,
		   const unsigned int n,
		   const RealVector& p_nTable ) const
{
    const Real thetaStep( M_PI / n );

    pTable.push_back( 0.0 );

    Real p_prev( 0.0 );
    unsigned int i( 1 );
    while( true )
    {
	const Real theta( thetaStep * i );

	Real p( this->p_theta_table( theta, r, r0, t, p_nTable ) );

	if( p < 0.0 )
	{
	    printf("drawTheta: p<0 %g\n", p );
	    p = 0.0;
	}

	const Real value( ( p_prev + p ) * 0.5 );
	pTable.push_back( *( pTable.end() - 1 ) + value );

//	printf("p %g %g %g\n", theta, pTable[i], p );

	if( //value < pTable[i] * std::numeric_limits<Real>::epsilon() ||
	    i >= n - 1 )
	{
	    break;   // pTable is valid in [0,i].
	}

	p_prev = p;
	++i;
    }

}


const Real 
FirstPassagePairGreensFunction::ip_theta( const Real theta,
					  const Real r, 
					  const Real r0, 
					  const Real t ) const
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in ( sigma, a )
	THROW_UNLESS( std::invalid_argument, r >= sigma && r < a );
	THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 || theta == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makep_nTable( p_nTable, r, r0, t );

    const Real p( ip_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}


const Real 
FirstPassagePairGreensFunction::idp_theta( const Real theta,
                                           const Real r, 
                                           const Real r0, 
                                           const Real t ) const
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in [ sigma, a ]
	THROW_UNLESS( std::invalid_argument, r >= sigma && r <= a );
	THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 || theta == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makedp_n_at_aTable( p_nTable, r0, t );

    const Real p( ip_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}

const Real 
FirstPassagePairGreensFunction::
ip_theta_n( const unsigned int n,
	    const RealVector& p_nTable, 
	    const RealVector& lgndTable1 ) const
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1( lgndTable1[n] );   // n-1
    const Real lgnd_n_p1( lgndTable1[n+2] ); // n+1
    
    return p_nTable[n] * ( lgnd_n_m1 - lgnd_n_p1 );// / ( 1.0 + 2.0 * n );
}


const Real 
FirstPassagePairGreensFunction::
ip_theta_table( const Real theta,
		const Real r, 
		const Real r0, 
		const Real t,	 
		const RealVector& p_nTable ) const
{
    const unsigned int tableSize( p_nTable.size() );

    const Real cos_theta( cos( theta ) );

    // LgndTable is offset by 1 to incorporate the n=-1 case.
    // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

    RealVector lgndTable( tableSize + 2 );
    lgndTable[0] = 1.0;  // n = -1
    gsl_sf_legendre_Pl_array( tableSize, cos_theta, &lgndTable[1] );

//    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
    const Real p( funcSum_all( boost::bind( &FirstPassagePairGreensFunction::
                                             ip_theta_n,
                                             this,
                                             _1, p_nTable, lgndTable ),
                                tableSize - 1 ) );

    return p;
}

const Real
FirstPassagePairGreensFunction::ip_theta_F( const Real theta,
					    const ip_theta_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r( params->r );
    const Real r0( params->r0 );
    const Real t( params->t );
    const RealVector& p_nTable( params->p_nTable );
    const Real value( params->value );

    return gf->ip_theta_table( theta, r, r0, t, p_nTable ) - value;
}



const Real 
FirstPassagePairGreensFunction::drawTheta( const Real rnd,
					   const Real r, 
					   const Real r0, 
					   const Real t ) const
{
    Real theta;

    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );
    THROW_UNLESS( std::invalid_argument, r >= sigma && r <= a );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 means no move.
    if( t == 0.0 )
    {
	return 0.0;
    }

    //const Real Dt6_r0( sqrt( 6.0 * getD() * t ) / r0 );
    //const Real high( 5 * Dt6_r0 < 1.0 ? asin( 5 * Dt6_r0 ) : M_PI );
    //const Real mean( Dt6_r0 < 1.0 ? asin( Dt6_r0 ) : M_PI );
    //printf("theta high %g %g %g\n", high,mean, t );
    const Real high( M_PI );

    RealVector p_nTable;

    if( r == geta() || r < 0.0 )
    {
	//puts("dp");
	makedp_n_at_aTable( p_nTable, r0, t );
    }
    else
    {
	makep_nTable( p_nTable, r, r0, t );
    }

    const Real ip_theta_pi( ip_theta_table( high, r, r0, t, p_nTable ) );

    ip_theta_params params = { this, r, r0, t, p_nTable, rnd * ip_theta_pi };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_theta_F ),
	    &params 
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, high );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	const int status( gsl_root_test_interval( low, high, 1e-15,
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawTheta: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    //printf("%d\n", i );

    theta = gsl_root_fsolver_root( solver );
    gsl_root_fsolver_free( solver );
    
    return theta;
}


//
// debug
//

const std::string FirstPassagePairGreensFunction::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getSigma() <<
	", a = " << this->geta() <<
	", kf = " << this->getkf() <<
	", h = " << this->geth() << std::endl;
    return ss.str();
}    
