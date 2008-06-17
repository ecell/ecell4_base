#if !defined( __FIRSTPASSAGENOCOLLISIONPAIRGREENSFUNCTION )
#define __FIRSTPASSAGENOCOLLISIONPAIRGREENSFUNCTION 

#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>

#include <gsl/gsl_roots.h>

#include "PairGreensFunction.hpp"

class FirstPassageNoCollisionPairGreensFunction
    :
    public PairGreensFunction
{

    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    static const Real MIN_T = 1e-18;

    static const unsigned int MAX_ORDER = 40;
    static const unsigned int MAX_ALPHA_SEQ = 1005;


public:
    
    FirstPassageNoCollisionPairGreensFunction( const Real D ); 
    
    virtual ~FirstPassageNoCollisionPairGreensFunction();

    const Real geta() const
    {
	return this->a;
    }

    void seta( const Real a );
    
    const Real drawTime( const Real rnd, const Real r0 ) const;

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
    
    const Real p_survival( const Real t,
			   const Real r0 ) const;

    const Real dp_survival( const Real t,
			    const Real r0 ) const;

    const Real p_int_r( const Real r,
			const Real t,
			const Real r0 ) const;

    const Real p_theta( const Real theta,
			const Real r, 
			const Real r0, 
			const Real t ) const;

    const Real ip_theta( const Real theta,
			 const Real r, 
			 const Real r0, 
			 const Real t ) const;

    const Real dp_theta( const Real theta,
                         const Real r, 
                         const Real r0, 
                         const Real t ) const;

    const Real idp_theta( const Real theta,
                          const Real r, 
                          const Real r0, 
                          const Real t ) const;


    const Real p_n( const Integer n, const Real r, 
		    const Real r0, const Real t ) const;

    const Real dp_n( const Integer n, const Real r0, const Real t ) const;


    const Real p_n_alpha( const unsigned int i,
			  const unsigned int n,
			  const Real r, 
			  const Real r0,
			  const Real t ) const;

    const Real dp_n_alpha( const unsigned int i,
                           const unsigned int n,
                           const Real r0,
                           const Real t ) const;
    
//    const unsigned int guess_maxi( const Real t ) const;

    // methods below are kept public for debugging purpose.

    const std::string dump() const;

protected:


    const Real p_theta_table( const Real theta,
			      const Real r, 
			      const Real r0, 
			      const Real t, 
			      const RealVector& p_nTable ) const;

    const Real ip_theta_table( const Real theta,
			       const Real r, 
			       const Real r0, 
			       const Real t,
			       const RealVector& p_nTable ) const;

    void makep_nTable( RealVector& p_nTable,
		       const Real r, 
		       const Real r0, 
		       const Real t ) const;
    
    void makedp_nTable( RealVector& p_nTable,
                        const Real r0, 
                        const Real t ) const;

    const Real p_theta_i( const unsigned int n,
			  const RealVector& p_nTable, 
			  const RealVector& lgndTable ) const;

    const Real ip_theta_i( const unsigned int n,
			   const RealVector& p_nTable, 
			   const RealVector& lgndTable1 ) const;

    
    struct p_survival_params
    { 
	const FirstPassageNoCollisionPairGreensFunction* const gf;
	const Real r0;
	const Real rnd;
    };

    static const Real 
    p_survival_F( const Real t,
		  const p_survival_params* const params );


    struct p_int_r_params
    { 
	const FirstPassageNoCollisionPairGreensFunction* const gf;
	const Real t;
	const Real r0;
//	const RealVector& num_r0Table;
	const Real rnd;
    };

    static const Real 
    p_int_r_F( const Real r,
	       const p_int_r_params* const params );

    struct ip_theta_params
    { 
	const FirstPassageNoCollisionPairGreensFunction* const gf;
	const Real r;
	const Real r0;
	const Real t;
	const RealVector& p_nTable;
	const Real value;
    };

    static const Real 
    ip_theta_F( const Real theta,
		const ip_theta_params* const params );


private:
    
    mutable boost::array<Integer,MAX_ORDER+1> alphaOffsetTable;
    mutable boost::array<RealVector,MAX_ORDER+1> alphaTable;
    //mutable std::vector<RealVector> alphaTable;

    Real a;

};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION 
