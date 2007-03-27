#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION 

#include <boost/tuple/tuple.hpp>

#include <gsl/gsl_roots.h>

#include "PairGreensFunction.hpp"



enum EventType
{
    REACTION = 0,
    ESCAPE = 1
};

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


    // true = reaction, false = shell.
    const EventType drawEventType( const Real rnd, 
				   const Real r0, 
				   const Real t ) const;

    
    const Real drawR( const Real rnd, 
		      const Real r0, 
		      const Real t ) const;
    
    const Real drawTheta( const Real rnd,
			  const Real r, 
			  const Real r0, 
			  const Real t ) const;
    
    
    const Real f_alpha0( const Real alpha ) const;
    const Real f_alpha0_aux( const Real alpha ) const;
  
    const Real f_alpha( const Real alpha, const Int n ) const;
    const Real f_alpha_aux( const Real alpha, const Int n ) const;

    
    const Real p_0( const Real t,
		    const Real r,
		    const Real r0 ) const;
    
    const Real p_survival( const Real t,
			   const Real r0 ) const;

    const Real p_leaves( const Real t,
			 const Real r0 ) const;

    const Real p_leavea( const Real t,
			 const Real r0 ) const;

    const Real p_int_r( const Real r,
			const Real t,
			const Real r0,
			const RealVector& num_r0Table ) const;

    const Real p_n( const Int n, const Real r, 
		    const Real r0, const Real t ) const;
    const Real p_n_alpha( const Real alpha,
			  const Int n, 
			  const Real r, 
			  const Real r0,
			  const Real t ) const;


protected:


    RealVector& getAlphaTable( const RealVector::size_type n ) const
    {
	if( this->alphaTable.size() <= n )
	{
	    this->alphaTable.resize( n+1 );
	}

	return this->alphaTable[n];
    }


    const Real alpha0_i( const Int i ) const;

    const Real alpha_i( const Int i, const Int n, 
			gsl_root_fsolver* const solver ) const;

    const Real p_0_i( const Real alpha,
		      const Real r,
		      const Real r0 ) const;

    const Real p_survival_i( const Real alpha,
			     const Real r0 ) const;

    const Real p_leavea_i( const Real alpha,
			   const Real r0 ) const;

    const Real p_leaves_i( const Real alpha,
			   const Real r0 ) const;

    const Real p_int_r_i( const Real r,
			  const Real alpha,
			  const Real r0,
			  const Real num_r0 ) const;

    const Real p_theta( const Real theta,
			const Real r, 
			const Real r0, 
			const Real t, 
			const RealVector& p_nTable ) const;
	
    void updateAlphaTable0( const Real t ) const;
    void updateAlphaTable( const Int n, const Real t ) const;
    void updateExpTable( const Real t ) const;
    void updatePsurvTable( const Real r0 ) const;
    void updateNum_r0Table( RealVector& num_r0Table,
			    const Real r0 ) const;

    void makep_nTable( const Real r, 
		      const Real r0, 
		      const Real t,
		      RealVector& p_nTable ) const;

    struct f_alpha0_aux_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Real value;
    };

    static const Real 
    f_alpha0_aux_F( const Real alpha,
		    const f_alpha0_aux_params* const params );


    struct f_alpha_aux_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Int n;
	Real value;
    };

    static const Real 
    f_alpha_aux_F( const Real alpha,
		   const f_alpha_aux_params* const params );

    
    struct p_survival_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Real r0;
	const Real rnd;
    };

    static const Real 
    p_survival_F( const Real t,
		  const p_survival_params* const params );


    struct p_int_r_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Real t;
	const Real r0;
	const Real psurv;
	const RealVector& num_r0Table;
	const Real rnd;
    };

    static const Real 
    p_int_r_F( const Real r,
	       const p_int_r_params* const params );
    
    const Real num_r0( const Real alpha,
		       const Real r0 ) const;

    static const Real P( const Int n, const Real x );
    static const Real Q( const Int n, const Real x );
    static const boost::tuple<Real,Real> P2( const Int n, const Real x );
    static const boost::tuple<Real,Real> Q2( const Int n, const Real x );



private:
    
    const Real h;
    const Real hsigma_p_1;

    mutable std::vector<RealVector> alphaTable;
    mutable RealVector expTable;
    mutable RealVector psurvTable;

    Real a;
    
    static const Real ALPHA_TOLERANCE = 1e-6;
    static const Real CUTOFF = 1e-8;

};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION 
