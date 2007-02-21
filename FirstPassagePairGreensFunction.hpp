#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION 

#include <iostream>
#include <math.h>

#include <gsl/gsl_integration.h>

#include "PairGreensFunction.hpp"

class FirstPassagePairGreensFunction
    :
    public PairGreensFunction
{
    
    
public:
    
    FirstPassagePairGreensFunction( const Real D, 
				    const Real kf, 
				    const Real Sigma );
    
    virtual ~FirstPassagePairGreensFunction();
    
    
    virtual const Real drawTime( const Real rnd, const Real r0,
				 const Real maxt ) const;
    
    virtual const Real drawR( const Real rnd, 
			      const Real r0, 
			      const Real t ) const;
    
    virtual const Real drawTheta( const Real rnd,
				  const Real r, 
				  const Real r0, 
				  const Real t ) const;
    
    
    const Real geth() const
    {
	return this->h;
    }


    const Real f_alpha_survival( const Real alpha, const Real a ) const;
    const Real f_alpha_survival_aux( const Real alpha, const Real a ) const;
    const Real f_alpha_survival_aux_df( const Real alpha, const Real a ) const;

    const Real alpha_survival_n( const Real a,
				 const Int n,
				 const Real lower,
				 const Real upper ) const;

  
    const Real f_alpha( const Real x, const Real a, const Int n ) const;

    const Real p_survival_i( const Real r,
			     const Real t,
			     const Real r0,
			     const Real alpha,
			     const Real a ) const;

private:
    
    
    struct f_alpha_survival_aux_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Real a;
	const Real value;
    };

    static const Real 
    f_alpha_survival_aux_F( const Real alpha,
			    const f_alpha_survival_aux_params* const params );
    static const Real 
    f_alpha_survival_aux_df_F( const Real alpha,
			       const f_alpha_survival_aux_params* const 
			       params );
    static void
    f_alpha_survival_aux_fdf_F( const Real alpha,
				const f_alpha_survival_aux_params* const 
				params,
				Real* const f, Real* const df );


    
    
private:
    
    const Real h;
//    const Real kD;
//    const Real alpha;
    
    static const Real P_CUTOFF = 1e-6;
    
};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION 
