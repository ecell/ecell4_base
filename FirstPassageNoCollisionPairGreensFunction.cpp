//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <algorithm>

#include <boost/bind.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_lambert.h>

#include "funcSum.hpp"
#include "freeFunctions.hpp"

#include "SphericalBesselGenerator.hpp"

#include "FirstPassageNoCollisionPairGreensFunction.hpp"



FirstPassageNoCollisionPairGreensFunction::
FirstPassageNoCollisionPairGreensFunction( const Real D ) 
    :
    PairGreensFunction( D, 0, 0 ),
    a( INFINITY )
{
    ; // do nothing
}

FirstPassageNoCollisionPairGreensFunction::
~FirstPassageNoCollisionPairGreensFunction()
{
    ; // do nothing
}

void FirstPassageNoCollisionPairGreensFunction::seta( const Real a )
{
    THROW_UNLESS( std::invalid_argument, a >= 0.0 );

    this->a = a;
}


const Real
FirstPassageNoCollisionPairGreensFunction::p_survival( const Real t,
                                                       const Real r0 ) const
{
    const Real D( getD() );
    const Real a( geta() );

    return p_survival_nocollision( t, r0, D, a );
}


const Real
FirstPassageNoCollisionPairGreensFunction::dp_survival( const Real t,
                                                        const Real r0 ) const
{
    const Real D( getD() );
    const Real a( geta() );

    return dp_survival_nocollision( t, r0, D, a );
}




const Real
FirstPassageNoCollisionPairGreensFunction::p_int_r( const Real r,
                                                    const Real t,
                                                    const Real r0 ) const
{
    const Real D( getD() );
    const Real a( geta() );

    const Real Dt( D * t );
    const Real asq( a * a );
    const Real a_r( 1.0 / a );
    const Real asq_r( a_r * a_r );

    const Real PIr0( M_PI * r0 );
    const Real PIr( M_PI * r );

    const Real r0_angle_factor( PIr0 * a_r );
    const Real r_angle_factor( PIr * a_r );
    const Real exp_factor( - Dt * M_PI * M_PI * asq_r );

    const unsigned int i_max( 
        std::max( static_cast<unsigned int>( 
                      ceil( sqrt( Dt * M_PI * M_PI + asq * 
                                  log( 1.0 / this->TOLERANCE ) / Dt ) * 
                            M_1_PI ) ), 
                  2u ) );

    Real p( 0.0 );
    unsigned int i( 1 );
    while( true )
    {
        Real sin_r;
        Real cos_r;
        sincos( r_angle_factor * i, &sin_r, &cos_r );

        const Real isq( i * i );

        const Real term1( exp( exp_factor * isq ) * 
                          sin( r0_angle_factor * i ) );
        const Real term2( a * sin_r - PIr * i * cos_r );
        const Real term( term1 * term2 / isq );
        
        p += term;

        if( i >= i_max )
        {
            break;
        }

        ++i;
    }

    const Real factor( M_2_PI / PIr0 );

    return p * factor;
}

const Real
FirstPassageNoCollisionPairGreensFunction::
p_survival_F( const Real t,
              const p_survival_params* params )
{
    const FirstPassageNoCollisionPairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    return rnd - gf->p_survival( t, r0 );
}


const Real
FirstPassageNoCollisionPairGreensFunction::
p_int_r_F( const Real r,
           const p_int_r_params* params )
{
    const FirstPassageNoCollisionPairGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    return gf->p_int_r( r, t, r0 ) - rnd;
}

/*
const unsigned int
FirstPassageNoCollisionPairGreensFunction::guess_maxi( const Real t ) const
{
    const Real D( getD() );
    const Real a( geta() );

    const Real Dt( D * t );
    const Real tolsq( this->TOLERANCE * this->TOLERANCE );
    const Real max_alpha( 1 / sqrt( gsl_sf_lambert_W0( 2 * Dt / tolsq ) 
                                    * tolsq  ) );
    
    return static_cast<unsigned int>( max_alpha * a / M_PI ) + 1;
}
*/


const Real 
FirstPassageNoCollisionPairGreensFunction::p_n_alpha( const unsigned int i,
						      const unsigned int n,
						      const Real r,
						      const Real r0, 
						      const Real t ) const
{
    const Real a( this->geta() );

    const Real mDt( - this->getD() * t );

    // j = a alpha -> alpha = j / a
    const Real aalpha( gsl_sf_bessel_zero_Jnu( static_cast<Real>( n ) + 0.5, 
                                               i + 1 ) );
    const Real alpha( aalpha / a );

    const Real term1( exp( mDt * alpha * alpha ) );

    const SphericalBesselGenerator& s( SphericalBesselGenerator::instance() );

    const Real jr(  s.j( n,   r * alpha ) );
    const Real jr0( s.j( n,   r0 * alpha ) );
    const Real ja2( s.j( n+1,   aalpha ) );

    const Real num( jr * jr0 );
    const Real den( ja2 * ja2 );

    const Real result( term1 * num / den );

    return result;
}


const Real 
FirstPassageNoCollisionPairGreensFunction::p_n( const Integer n,
                                                const Real r,
                                                const Real r0, 
                                                const Real t ) const
{
    const Real p( funcSum( 
                      boost::bind( &FirstPassageNoCollisionPairGreensFunction::
                                   p_n_alpha,
                                   this,
                                   _1, n, r, r0, t ),
                      this->MAX_ALPHA_SEQ ) );

    return p;
}

void
FirstPassageNoCollisionPairGreensFunction::makep_nTable( RealVector& p_nTable,
                                                         const Real r, 
                                                         const Real r0, 
                                                         const Real t ) const
{
    const Real a( geta() );

    p_nTable.clear();

    const Real factor( 1.0 / ( 2.0 * M_PI * gsl_pow_3( a ) ) ); 

    const Real p_0( this->p_n( 0, r, r0, t ) * factor );
    p_nTable.push_back( p_0 );
    //printf("0 p_n %18.18g\n", p_0 );
    const Real threshold( fabs( p_0 * THETA_TOLERANCE * 1e-1  ) );

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
	//printf("%d p_n %18.18g\n", n, p_n );

	p_nTable.push_back( p_n );

	const Real p_n_abs( fabs( p_n ) );
	// truncate when converged enough.
	if( p_n_abs <= threshold &&
            p_n_prev_abs <= threshold &&
	    p_n_abs <= p_n_prev_abs )
	{
	    break;
        }
	
	++n;

	if( n >= this->MAX_ORDER )
	{
	    std::cerr << "p_n didn't converge." << std::endl;
	    break;
	}
	
	p_n_prev_abs = p_n_abs;
    }

}



const Real 
FirstPassageNoCollisionPairGreensFunction::
p_theta_i( const unsigned int n,
	   const RealVector& p_nTable, const RealVector& lgndTable ) const
{
    return p_nTable[n] * lgndTable[n] * ( 2 * n + 1 );
}

const Real
FirstPassageNoCollisionPairGreensFunction::
p_theta_table( const Real theta,
	       const Real r, 
	       const Real r0, 
	       const Real t, 
	       const RealVector& p_nTable ) const
{
    const unsigned int tableSize( p_nTable.size() );

    RealVector lgndTable( tableSize );

    Real sin_theta;
    Real cos_theta;
    sincos( theta, &sin_theta, &cos_theta );
    gsl_sf_legendre_Pl_array( tableSize-1, cos_theta, &lgndTable[0] );

    const Real p( funcSum_all( 
                      boost::bind( &FirstPassageNoCollisionPairGreensFunction::
                                   p_theta_i,
                                   this,
                                   _1, p_nTable, lgndTable ),
                      tableSize ) );
    
    return p * sin_theta;
}


const Real 
FirstPassageNoCollisionPairGreensFunction::p_theta( const Real theta,
                                                    const Real r, 
                                                    const Real r0, 
                                                    const Real t ) const 
{
    {
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
	THROW_UNLESS( std::invalid_argument, r >= 0 && r < a );
	THROW_UNLESS( std::invalid_argument, r0 >= 0 && r0 < a );
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
FirstPassageNoCollisionPairGreensFunction::ip_theta( const Real theta,
                                                     const Real r, 
                                                     const Real r0, 
                                                     const Real t ) const
{
    {
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in ( sigma, a )
	THROW_UNLESS( std::invalid_argument, r >= 0.0 && r < a );
	THROW_UNLESS( std::invalid_argument, r0 >= 0.0 && r0 < a );
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
FirstPassageNoCollisionPairGreensFunction::
ip_theta_i( const unsigned int n,
	    const RealVector& p_nTable, 
	    const RealVector& lgndTable1 ) const
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1( lgndTable1[n] );   // n-1
    const Real lgnd_n_p1( lgndTable1[n+2] ); // n+1
    
    return p_nTable[n] * ( lgnd_n_m1 - lgnd_n_p1 );// / ( 1.0 + 2.0 * n );
}


const Real 
FirstPassageNoCollisionPairGreensFunction::
ip_theta_table( const Real theta,
		const Real r, 
		const Real r0, 
		const Real t,	 
		const RealVector& p_nTable ) const
{
    const unsigned int tableSize( p_nTable.size() );

    RealVector pTable;
    pTable.reserve( tableSize );

    const Real cos_theta( cos( theta ) );

    // LgndTable is offset by 1 to incorporate the n=-1 case.
    // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

    RealVector lgndTable1( tableSize + 2 );
    lgndTable1[0] = 1.0;  // n = -1
    gsl_sf_legendre_Pl_array( tableSize, cos_theta, &lgndTable1[1] );


    const Real p( funcSum_all( boost::bind( 
                               &FirstPassageNoCollisionPairGreensFunction::
                               ip_theta_i,
                               this,
                               _1, p_nTable, lgndTable1 ),
			   tableSize ) );

    return p;
}



const Real
FirstPassageNoCollisionPairGreensFunction::
ip_theta_F( const Real theta,
            const ip_theta_params* params )
{
    const FirstPassageNoCollisionPairGreensFunction* const gf( params->gf ); 
    const Real r( params->r );
    const Real r0( params->r0 );
    const Real t( params->t );
    const RealVector& p_nTable( params->p_nTable );
    const Real value( params->value );

    return gf->ip_theta_table( theta, r, r0, t, p_nTable ) - value;
}


const Real 
FirstPassageNoCollisionPairGreensFunction::dp_n_alpha( const unsigned int i,
                                                       const unsigned int n,
                                                       const Real r0, 
                                                       const Real t ) const
{
    const Real a( this->geta() );

    const Real mDt( - this->getD() * t );

    const Real 
        aalpha( gsl_sf_bessel_zero_Jnu( static_cast<Real>( n ) + 0.5, i + 1 ) );
    const Real alpha( aalpha / a );

    const Real term1( exp( mDt * alpha * alpha ) * alpha );

    const SphericalBesselGenerator& s( SphericalBesselGenerator::instance() );

    const Real jr0( s.j( n,   r0 * alpha ) );
    const Real ja2( s.j( n+1,   aalpha ) );

    const Real result( term1 * jr0 / ja2 );

    return result;
}



const Real 
FirstPassageNoCollisionPairGreensFunction::dp_n( const Integer n,
                                                 const Real r0, 
                                                 const Real t ) const
{
    const Real 
        p( funcSum( boost::bind( &FirstPassageNoCollisionPairGreensFunction::
                                 dp_n_alpha,
                                 this,
                                 _1, n, r0, t ),
                    this->MAX_ALPHA_SEQ ) );

    return p;
}


void
FirstPassageNoCollisionPairGreensFunction::
makedp_nTable( RealVector& p_nTable,
               const Real r0, 
               const Real t ) const
{
    p_nTable.clear();

    const Real factor( - getD() / ( 2.0 * M_PI * gsl_pow_3( a ) ) );

    const Real p_0( this->dp_n( 0, r0, t ) * factor );
    p_nTable.push_back( p_0 );

    const Real threshold( fabs( THETA_TOLERANCE * p_0 * 1e-1 ) );

    Real p_n_prev_abs( fabs( p_0 ) );
    unsigned int n( 1 );
    while( true )
    {
	Real p_n( this->dp_n( n, r0, t ) * factor );

	if( ! std::isfinite( p_n ) )
	{
	    std::cerr << "makedp_nTable: invalid value; " <<
		p_n << "( n= " << n << ")." << std::endl;
//	    p_n = 0.0;
	    break;
	}
	//printf("dp_n %g\n",p_n );

	p_nTable.push_back( p_n );

	const Real p_n_abs( fabs( p_n ) );
	// truncate when converged enough.
	if( p_n_abs <= threshold &&
            p_n_prev_abs <= threshold &&
	    p_n_abs <= p_n_prev_abs )
	{
	    break;
	}
	
	++n;

	if( n >= this->MAX_ORDER )
	{
	    //std::cerr << "dp_n didn't converge." << std::endl;
	    break;
	}
	
	p_n_prev_abs = p_n_abs;
    }

}

const Real 
FirstPassageNoCollisionPairGreensFunction::dp_theta( const Real theta,
                                                     const Real r, 
                                                     const Real r0, 
                                                     const Real t ) const 
{
    {
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );

        // r \in [ sigma, a ]  ;  unlike p_theta,
        // defined at r == sigma and r == a.
	THROW_UNLESS( std::invalid_argument, r >= 0.0 && r <= a );
	THROW_UNLESS( std::invalid_argument, r0 >= 0.0 && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makedp_nTable( p_nTable, r0, t );

    const Real p( p_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}

const Real 
FirstPassageNoCollisionPairGreensFunction::idp_theta( const Real theta,
                                                      const Real r, 
                                                      const Real r0, 
                                                      const Real t ) const
{
    {
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in [ sigma, a ]
	THROW_UNLESS( std::invalid_argument, r >= 0.0 && r <= a );
	THROW_UNLESS( std::invalid_argument, r0 >= 0.0 && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 || theta == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makedp_nTable( p_nTable, r0, t );

    const Real p( ip_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}

const Real 
FirstPassageNoCollisionPairGreensFunction::drawTime( const Real rnd, 
                                                     const Real r0 ) const
{
   const Real a( this->geta() );

   THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
   THROW_UNLESS( std::invalid_argument, r0 >= 0.0 && r0 <= a );

   if( r0 == a || a == 0.0 )
   {
       return 0.0;
   }

   Real low( 1e-6 );
   Real high( 1.0 );

   p_survival_params params = { this, r0, rnd };

   gsl_function F = 
       {
           reinterpret_cast<typeof(F.function)>( &p_survival_F ),
           &params 
       };

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    while( GSL_FN_EVAL( &F, high ) < 0.0 )
    {
	high *= 10;
	printf( "drawTime: adjusting high: %g\n", high );
	if( fabs( high ) >= 1e10 )
	{
	    std::cerr << "Couldn't adjust high. F(" << high <<
		") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
		", " << dump() << std::endl;
	    throw std::exception();
	}
    }

    Real low_value( GSL_FN_EVAL( &F, low ) );
    while( low_value > 0.0 )
    {
	low *= .1;

        const Real low_value_new( GSL_FN_EVAL( &F, low ) );

	printf( "drawTime: adjusting low: %g, F = %g\n", low, low_value_new );

	if( fabs( low ) <= this->MIN_T || 
            fabs( low_value - low_value_new ) < TOLERANCE ) 
	{
	    std::cerr << "Couldn't adjust low.  Returning MIN_T (= "
		      << this->MIN_T << "); F(" << low <<
		") = " << GSL_FN_EVAL( &F, low ) << "; r0 = " << r0 << ", "
		      << dump() << std::endl;
	    return this->MIN_T;
	}

        low_value = low_value_new;
    }

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

	const int status( gsl_root_test_interval( low, high, this->MIN_T, 
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawTime: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    // printf("%d\n", i );


    Real t( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return t;
}



const Real 
FirstPassageNoCollisionPairGreensFunction::drawR( const Real rnd, 
                                                  const Real r0, 
                                                  const Real t ) const
{
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= 0.0 && r0 < a );

    if( t == 0.0 )
    {
	return r0;
    }

    const Real psurv( p_survival( t, r0 ) );

    p_int_r_params params = { this, t, r0, rnd * psurv };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};

    Real low( 0.0 );
    Real high( a );

//    const Real lowvalue( GSL_FN_EVAL( &F, low  ) );
    const Real highvalue( GSL_FN_EVAL( &F, high ) );

    // No initial range guess, except the negative value check below,
    // as evaluation of p_int_r in this GF seems pretty robust.

    if( highvalue < 0.0 )
    {
	printf( "drawR: highvalue < 0.0 (%g). returning a.\n", highvalue );
	return a;
    }


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
    
const Real 
FirstPassageNoCollisionPairGreensFunction::drawTheta( const Real rnd,
                                                      const Real r, 
                                                      const Real r0, 
                                                      const Real t ) const
{
    Real theta;

    const Real a( this->geta() );

    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= 0.0 && r0 < a );
    THROW_UNLESS( std::invalid_argument, r >= 0.0 && r <= a );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 means no move.
    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    if( r == geta() || r < 0.0 )
    {
	//puts("dp");
	makedp_nTable( p_nTable, r0, t );
    }
    else
    {
	makep_nTable( p_nTable, r, r0, t );
    }

    // root finding with the integrand form.

    const Real ip_theta_pi( ip_theta_table( M_PI, r, r0, t, p_nTable ) );

    ip_theta_params params = { this, r, r0, t, p_nTable, rnd * ip_theta_pi };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_theta_F ),
	    &params 
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, M_PI );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	const int status( gsl_root_test_interval( low, high, 1e-11,
						  THETA_TOLERANCE ) );

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

const std::string FirstPassageNoCollisionPairGreensFunction::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() <<
	", a = " << this->geta() << std::endl;
    return ss.str();
}    
