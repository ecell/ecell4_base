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

  virtual const Real drawTime( const Real rnd, const Real r0,
			       const Real maxt ) const;

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


  

private:

  struct p_reaction_params 
  { 
    const PlainPairGreensFunction* const gf;
    const Real r0;
    const Real rnd;
  };

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

  static const Real p_reaction_F( const Real tsqrt, 
				  const p_reaction_params* params );
  static const Real p_reaction_deriv_F( const Real tsqrt,
					const p_reaction_params* params );
  static void p_reaction_fdf_F( const Real tsqrt,
				const p_reaction_params* params,
				Real* const f, Real* const df );

  const Real p_reaction( const Real tsqrt, const Real r0 ) const;
  const Real p_reaction_deriv( const Real tsqrt, const Real r0 ) const;
  void p_reaction_fdf( const Real tsqrt, const Real r0,
		       Real* const f, Real* const df ) const;

  const Real p_irr_radial( const Real r, const Real r0, const Real t ) const;

  const Real p_free( const Real r, const Real r0, 
		     const Real theta, const Real t ) const;

  const Real p_corr( const Real r, const Real r0, 
		     const Real theta, const Real t ) const;

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
      gsl_integration_workspace* const workspace, const Real err ) const;


private:

  const Real kD;
  const Real alpha;

  static const Real P_CUTOFF = 1e-6;
  static const Real H = 3.0;

};



#endif // __PLAINPAIRGREENSFUNCTION 
