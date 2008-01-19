#if !defined( __PLAINPAIRGREENSFUNCTION )
#define __PLAINPAIRGREENSFUNCTION 

#include <iostream>
#include <cmath>

#include <gsl/gsl_integration.h>

#include "PairGreensFunction.hpp"

class PlainPairGreensFunction
    :
    public PairGreensFunction
{
    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    static const Real MIN_T = 1e-12;

    static const unsigned int MAX_ORDER = 45;

    
public:
    
    PlainPairGreensFunction( const Real D, const Real kf, const Real Sigma );
    
    ~PlainPairGreensFunction();
    
    
    const Real drawTime( const Real rnd, const Real r0 ) const;
    
    const Real drawR( const Real rnd, 
		      const Real r0, 
		      const Real t ) const;
    
    const Real drawTheta( const Real rnd,
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
    


    const Real p_tot( const Real r, const Real r0, 
		      const Real theta, const Real time ) const;
    
    const Real p_theta( const Real r, const Real r0, 
                        const Real theta, const Real time ) const;

    const Real ip_theta( const Real r, const Real r0, 
                         const Real theta, const Real time ) const;

    const Real p_free( const Real theta, const Real r, const Real r0, 
		       const Real t ) const;

    const Real ip_free( const Real theta, 
                        const Real r, 
                        const Real r0, 
                        const Real t ) const;
    
    const Real p_corr( const Real r, const Real r0, 
		       const Real theta, const Real t ) const;

    const std::string dump() const;

private:
    
    struct p_reaction_params 
    { 
	const PlainPairGreensFunction* const gf;
	const Real r0;
	const Real rnd;
    };

    struct p_theta_params 
    { 
	const PlainPairGreensFunction* const gf;
	const Real r;
	const Real r0;
	const Real t;
        const RealVector& RnTable;
	const Real value;
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

    static const Real ip_theta_F( const Real theta,
                                  const p_theta_params* params );

    const Real p_reaction( const Real t, const Real r0 ) const;

    static const Real p_corr_R( const Real u, 
				const p_corr_R_params* const params );
    static const Real p_corr_R2( const Real u, 
				 const p_corr_R2_params* const params );
    
    const Real p_corr_n( const unsigned int n, const RealVector& RnTable, 
                         const RealVector& lgndTable ) const;

    const Real ip_corr_n( const unsigned int n, const RealVector& RnTable, 
                          const RealVector& lgndTable ) const;

    const Real 
    p_corr_table( const Real theta, const Real r, const Real r0,
                   const Real t, const RealVector& RnTable ) const;

    const Real 
    ip_corr_table( const Real theta, const Real r, const Real r0,
                   const Real t, const RealVector& RnTable ) const;
    
    const Real p_theta_table( const Real r, const Real r0, 
                              const Real theta, const Real time,
                              const RealVector& RnTable ) const;

    const Real ip_theta_table( const Real r, const Real r0, 
                               const Real theta, const Real time,
                               const RealVector& RnTable ) const;

    const Real 
    p_corr_table( const Real theta, const Real r, const Real r0,
                  const Real t, const RealVector& RnTable );
    

    void makeRnTable( RealVector& RnTable,
                      const Real r,
                      const Real r0,
                      const Real t ) const;

    const Real 
    Rn( const Integer order, const Real r, const Real r0, const Real t,
	gsl_integration_workspace* const workspace, const Real err ) const;
    
    
private:
    
    const Real kD;
    const Real alpha;
    
    static const Real P_CUTOFF = 1e-6;
    static const Real H = 3.0;
    
};



#endif // __PLAINPAIRGREENSFUNCTION 
