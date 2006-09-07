#if !defined( __PLAINPAIRGREENSFUNCTION )
#define __PLAINPAIRGREENSFUNCTION 

#include <iostream>
#include <math.h>

#include <gsl/gsl_integration.h>

#include "PairGreensFunction.hpp"

class PlainPairGreensFunction
  :
  public PairGreensFunction
{


public:

  PlainPairGreensFunction( const Real D, const Real kf, const Real Sigma );

  virtual ~PlainPairGreensFunction();


  virtual const Real p_tot( const Real r, const Real r0, 
			    const Real theta, const Real time ) const;

  virtual const Real drawNextReactionTime( const Real rnd, const Real r0 );

  virtual const Real drawR( const Real rnd, 
			    const Real r0, 
			    const Real t ) const;

  virtual const Real drawTheta( const Real rnd,
				const Real r, 
				const Real r0, 
				const Real t ) const;


  const Real getkD() const
  {
    return this->kD;
  }

  const Real getalpha() const
  {
    return this->alpha;
  }


  
  const Real p_free(  const Real r, const Real r0, 
		      const Real theta, const Real t ) const;

  const Real p_corr( const Real r, const Real r0, 
		     const Real theta, const Real t ) const;

  const Real p_irr_radial( const Real r, const Real r0, 
			   const Real t ) const;

  const Real intt_p_irr_radial( const Real r, 
				const Real r0, 
				const Real t ) const;


private:

  struct p_corr_R_params 
  { 
    int order;
    const Real r;
    const Real r0; 
    const Real t; 
    const Real Sigma; 
    const Real D;
    const Real kf;
  };

  struct p_corr_R2_params 
  { 
    int order;
    int order_max;
    const Real r;
    const Real r0; 
    const Real theta;
    const Real t; 
    const Real Sigma; 
    const Real D;
    const Real kf;
  };

  const Real p_survival( const Real t, const Real r0 ) const;
  const Real p_survival_deriv( const Real t, const Real r0 ) const;


  static const Real p_corr_R( const Real u, 
			      const p_corr_R_params* const params );
  static const Real p_corr_R2( const Real u, 
			       const p_corr_R2_params* const params );


  const Real p_tot_RnTable( const Real r, const Real r0, 
			    const Real theta, const Real time,
			    const RealVector& RnTable ) const;

  static const Real 
  p_corr_RnTable( const Real theta, const Real r, const Real r0,
		  const Real t, const RealVector& RnTable );

  const Real 
  Rn( const Int order, const Real r, const Real r0, const Real t,
      gsl_integration_workspace* workspace ) const;

  const Real getMinR() const
  {
    return getSigma();   // need to consider r0?
  }

  const Real getMaxR( const Real t ) const
  {
    // [ Sigma, Sigma + 2 H sqrt( 6 D t ) ];
    // should actually be [ Sigma, r0(t) + H sqrt( 6 D t ) ]
    return getMinR() + (this->H*2) * sqrt( 6.0 * getD() * t );
  }



private:

  //  const Real kD;
  //  const Real alpha;
   Real kD;
   Real alpha;

  static const Real P_CUTOFF = 1e-6;
  static const Real H = 4.0;

};



#endif // __PLAINPAIRGREENSFUNCTION 
