#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION 

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

    const Real geth() const
    {
	return this->h;
    }

    const Real geta() const
    {
	return this->a;
    }

    void seta( Real a );
    
    
    const Real drawTime( const Real rnd, const Real r0 ) const;
    
    const Real drawR( const Real rnd, 
		      const Real r0, 
		      const Real t ) const;
    
    const Real drawTheta( const Real rnd,
			  const Real r, 
			  const Real r0, 
			  const Real t ) const;
    
    
    const Real f_alpha_survival( const Real alpha ) const;
    const Real f_alpha_survival_aux( const Real alpha ) const;
    const Real f_alpha_survival_aux_df( const Real alpha ) const;

    const Real alpha_survival_n( const Int n ) const;

  
    const Real f_alpha( const Real x, const Int n ) const;

    const Real p_survival( const Real t,
			   const Real r0 ) const;

    const Real p_leaves( const Real t,
			 const Real r0 ) const;

    const Real p_leavea( const Real t,
			 const Real r0 ) const;


    const Real p_survival_i( const Real alpha,
			     const Real r0 ) const;

    const Real p_leavea_i( const Real alpha,
			   const Real r0 ) const;

    const Real p_leaves_i( const Real alpha,
			   const Real r0 ) const;

    const Real asratio( const Real alpha,
			const Real r0 ) const;



    const RealVector& getAlphaTable() const
    {
	return this->alphaTable;
    }


protected:

    void updateAlphaTable( const Real t ) const;
    void updateExpTable( const Real t ) const;
    
    struct f_alpha_survival_aux_params
    { 
	const FirstPassagePairGreensFunction* const gf;
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

    struct p_survival_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Real r0;
	const Real rnd;
    };

    static const Real 
    p_survival_F( const Real t,
		  const p_survival_params* const params );
    
    
private:
    
    const Real h;
    const Real hsigma_p_1;

    mutable RealVector alphaTable;
    mutable RealVector expTable;

    Real a;
    
    static const Real CUTOFF = 1e-8;

    static const Real ALPHA_CUTOFF = 1e-8;
    
};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION 
