#if !defined( __PLAINPAIRGREENSFUNCTION )
#define __PLAINPAIRGREENSFUNCTION 

#include <iostream>
#include <cmath>

#include <gsl/gsl_integration.h>

#include "PairGreensFunction.hpp"

class BasicPairGreensFunction
    :
    public PairGreensFunction
{
    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    static const Real MIN_T = 1e-12;

    static const unsigned int MAX_ORDER = 65;

    static const Real H = 4.0;
    

    
public:
    
    BasicPairGreensFunction( const Real D, const Real kf, const Real Sigma );
    
    ~BasicPairGreensFunction();
    
    
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
    
    const Real p_reaction( const Real t, const Real r0 ) const;
    const Real p_survival( const Real t, const Real r0 ) const;
    const Real p_int_r( const Real r, 
                        const Real t, 
                        const Real r0 ) const;

//    const Real p_int_r_max( const Real t, const Real r0 ) const;

    
    const Real p_theta( const Real theta, const Real r, const Real r0, 
                        const Real time ) const;

    const Real ip_theta( const Real theta, const Real r, const Real r0, 
                         const Real time ) const;

    const Real p_free( const Real theta, const Real r, const Real r0, 
		       const Real t ) const;

    const Real ip_free( const Real theta, 
                        const Real r, 
                        const Real r0, 
                        const Real t ) const;
    
    const Real p_corr( const Real theta, const Real r, const Real r0, 
		       const Real t ) const;

    const Real ip_corr( const Real theta, const Real r, const Real r0, 
                        const Real t ) const;

    const std::string dump() const;

private:
    
    struct p_reaction_params 
    { 
	const BasicPairGreensFunction* const gf;
	const Real r0;
	const Real rnd;
    };

    struct p_int_r_params 
    { 
	const BasicPairGreensFunction* const gf;
	const Real t;
	const Real r0;
	const Real rnd;
    };

    struct p_theta_params 
    { 
	const BasicPairGreensFunction* const gf;
	const Real r;
	const Real r0;
	const Real t;
        const RealVector& RnTable;
	const Real value;
    };
    
    struct p_corr_R_params 
    { 
	const BasicPairGreensFunction* const gf;
	unsigned int n;
	const Real r;
	const Real r0; 
	const Real t; 
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

    static const Real p_int_r_F( const Real r,
                                 const p_int_r_params* const params );

    static const Real p_corr_R_F( const Real alpha, 
                                  const p_corr_R_params* const params );
//    static const Real p_corr_R2( const Real alpha, 
//				 const p_corr_R2_params* const params );


    const Real p_corr_R( const Real alpha,
                         const unsigned int n,
                         const Real r,
                         const Real r0,
                         const Real t ) const;

    
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
    Rn( const unsigned int order, const Real r, const Real r0, const Real t,
	gsl_integration_workspace* const workspace, const Real tol ) const;
    
    
private:
    
    const Real kD;
    const Real alpha;
    
};



#endif // __PLAINPAIRGREENSFUNCTION 
