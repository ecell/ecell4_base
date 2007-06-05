#include "FreePairGreensFunction.hpp"



const Real 
FreePairGreensFunction::drawR( const Real rnd, 
                               const Real r0, 
                               const Real t ) const
{

}

const Real 
FreePairGreensFunction::drawTheta( const Real rnd,
                                   const Real r, 
                                   const Real r0, 
                                   const Real t ) const
{

}


const Real 
FreePairGreensFunction::ip_r( const Real r, 
                              const Real r0, 
                              const Real t ) const
{

}
    
const Real 
FreePairGreensFunction::p_theta( const Real theta,
                                 const Real r,
                                 const Real r0,
                                 const Real t ) const
{
    const Real D( getD() );
    const Real Dt4( 4.0 * D * t );
    const Real Dt4PI( Dt4 * M_PI );

    const Real p( ( 1.0 / std::sqrt( gsl_pow_3( Dt4PI ) ) ) *
                  exp( - ( r * r + r0 * r0 - 2.0 * r * r0 * cos( theta ) ) 
                       / Dt4 ) );

    return p * sin( theta );
}

const Real 
FreePairGreensFunction::ip_theta( const Real theta,
                                  const Real r,
                                  const Real r0,
                                  const Real t ) const
{
    const Real Dt( getD() * t );
    const Real Dt2( Dt + Dt );
    const Real rr0( r * r0 );

    const Real rr0_over_2Dt( rr0 / ( Dt2 ) );

    const Real rsqr0sq_over_4Dt( ( r * r + r0 * r0 ) / ( Dt2 + Dt2 ) );

    const Real term1( expm1( rr0_over_2Dt 
                             - rsqr0sq_over_4Dt ) );
    const Real term2( expm1( rr0_over_2Dt * cos( theta ) 
                             - rsqr0sq_over_4Dt ) );

    const Real den( 4.0 * sqrt( M_PI * M_PI * M_PI * Dt ) * rr0 );

    return ( term1 - term2 ) / den;
}

/*
const Real 
FirstPassagePairGreensFunction::dp_theta_free( const Real theta,
                                               const Real r,
                                               const Real r0,
                                               const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real Dt4( 4.0 * Dt );
    const Real rsq( r * r );

    Real sin_theta;
    Real cos_theta;
    sincos( theta, &sin_theta, &cos_theta );
    const Real rr0costheta( r * r0 * cos_theta );

    const Real num1( r * D * 
                    exp( - ( rsq + r0 * r0 - 2.0 * rr0costheta ) / 
                         Dt4 ) );

    const Real num2( - rsq + Dt4 + rr0costheta );

    const Real den( 4.0 * sqrt( M_PI * gsl_pow_5( Dt ) ) );

    return ( num1 * num2 / den ) * sin_theta;
}

const Real 
FirstPassagePairGreensFunction::idp_theta_free( const Real theta,
                                                const Real r,
                                                const Real r0,
                                                const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real Dt2( Dt + Dt );
    const Real rr0( r * r0 );
    const Real rsq( r * r );
    const Real rr0_over_2Dt( rr0 / ( Dt2 ) );
    const Real rsqr0sq_over_4Dt( ( rsq + r0 * r0 ) / ( Dt2 + Dt2 ) );

    const Real cos_theta( cos( theta ) );

    const Real num1( exp( rr0_over_2Dt - rsqr0sq_over_4Dt ) *
                      ( - rsq + rr0 + Dt2 ) );
    const Real num2( exp( rr0_over_2Dt * cos_theta - rsqr0sq_over_4Dt ) *
                      ( - rsq + rr0 * cos_theta + Dt2 ) );

    const Real den( 2.0 * sqrt( M_PI * Dt * Dt * Dt ) * r0 );

    return D * ( num1 - num2 ) / den;
}
*/
