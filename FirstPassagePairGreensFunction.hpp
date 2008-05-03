#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION_HPP )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION_HPP 

#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>

#include <gsl/gsl_roots.h>

#include "PairGreensFunction.hpp"


class FirstPassagePairGreensFunction
    :
    public PairGreensFunction
{

    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    static const Real MIN_T = 1e-12;

    static const unsigned int MAX_ORDER = 100;
    static const unsigned int MAX_ALPHA_SEQ = 1000;


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

    void seta( const Real a );
    
    const Real drawTime( const Real rnd, const Real r0 ) const;

    const boost::tuple<Real,EventType> 
    drawTime2( const Real rnd1, const Real rnd2, const Real r0,
               const Real t_guess = 1e-3 ) const;


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
  
    const Real f_alpha( const Real alpha, const Integer n ) const;
    const Real f_alpha_aux( const Real alpha, const Integer n ) const;

    
    const Real p_0( const Real t,
		    const Real r,
		    const Real r0 ) const;
    
    const Real p_survival( const Real t,
			   const Real r0 ) const;

    const Real p_survival_table( const Real t,
				 const Real r0,
				 const RealVector& psurvTable ) const;

    const Real dp_survival( const Real t,
			    const Real r0 ) const;

    const Real leaves( const Real t,
		       const Real r0 ) const;

    const Real leavea( const Real t,
		       const Real r0 ) const;

    const Real p_leaves( const Real t,
			 const Real r0 ) const;

    const Real p_leavea( const Real t,
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

    const Real dp_n_at_a( const Integer n, const Real r0, const Real t ) const;


    const Real p_n_alpha( const unsigned int i,
			  const unsigned int n,
			  const Real r, 
			  const Real r0,
			  const Real t ) const;
//                          const RealVector& jas1Table ) const;

    const Real dp_n_alpha_at_a( const unsigned int i,
				const unsigned int n,
				const Real r0,
				const Real t ) const;

    // methods below are kept public for debugging purpose.

    const std::string dump() const;

    const unsigned int alphaOffset( const unsigned int n ) const;

    const Real alpha0_i( const Integer i ) const;

    const Real alpha_i( const Integer i, const Integer n, 
			gsl_root_fsolver* const solver ) const;


protected:

    void clearAlphaTable() const;


    RealVector& getAlphaTable( const size_t n ) const
    {
	return this->alphaTable[n];
    }

    const Real getAlpha( const size_t n, const RealVector::size_type i ) const
    {
	RealVector& alphaTable( this->alphaTable[n] );
        const RealVector::size_type oldSize( alphaTable.size() );

	if( oldSize <= i )
	{
	    alphaTable.resize( i+1 );
	    const unsigned int offset( alphaOffset( n ) );

            const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
            gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

            for( RealVector::size_type m( oldSize ); m <= i; ++m )
            {
                alphaTable[m] = alpha_i( m + offset, n, solver );
            }

            gsl_root_fsolver_free( solver );
	}

	return alphaTable[i];

    }

    const Real getAlpha0( const RealVector::size_type i ) const
    {
	RealVector& alphaTable( this->alphaTable[0] );
	
        const RealVector::size_type oldSize( alphaTable.size() );

	if( oldSize <= i )
	{
	    alphaTable.resize( i+1 );

            for( RealVector::size_type m( oldSize ); m <= i; ++m )
            {
                alphaTable[m] = alpha0_i( m );
            }
	}

	return alphaTable[i];
    }



    const Real p_0_i( const Real alpha,
		      const Real r,
		      const Real r0 ) const;

    const Real p_survival_i( const Real alpha,
			     const Real r0 ) const;

    const Real dp_survival_i( const Real alpha,
			      const Real r0 ) const;

    const Real leavea_i( const Real alpha,
			 const Real r0 ) const;

    const Real leaves_i( const Real alpha,
			 const Real r0 ) const;

    const Real p_leavea_i( const Real alpha,
			   const Real r0,
                           const Real pleave_factor ) const;

    const Real p_leaves_i( const Real alpha,
			   const Real r0,
                           const Real pleave_factor ) const;

    const Real p_survival_den( const Real alpha,
                               const Real r0 ) const;

    const Real p_int_r_i( const Real r,
			  const Real alpha,
			  const Real r0,
			  const Real num_r0 ) const;

    const Real p_int_r_table( const Real r,
			      const Real t,
			      const Real r0,
			      const RealVector& num_r0Table ) const;

    const Real ip_theta_table( const Real theta,
			       const Real r, 
			       const Real r0, 
			       const Real t,
			       const RealVector& p_nTable ) const;

    const Real dp_theta_at_a( const Real theta,
			      const Real r0, 
			      const Real t ) const;


    const Real p_theta_table( const Real theta,
			      const Real r, 
			      const Real r0, 
			      const Real t, 
			      const RealVector& p_nTable ) const;

    void make_p_thetaTable( RealVector& pTable,
			    const Real r, 
			    const Real r0, 
			    const Real t,
			    const unsigned int n,
			    const RealVector& p_nTable ) const;

    const Real p_0_i_exp( const unsigned int i,
			  const Real t,
			  const Real r,
			  const Real r0 ) const;

    const Real p_survival_i_exp( const unsigned int i,
				 const Real t,
				 const Real r0 ) const;

    const Real p_survival_i_exp_table( const unsigned int i,
				       const Real t,
				       const Real r0,
				       const RealVector& psurvTable ) const;

    const Real dp_survival_i_exp( const unsigned int i,
				  const Real alpha,
				  const Real r0 ) const;

    const Real leavea_i_exp( const unsigned int i,
			     const Real alpha,
			     const Real r0 ) const;

    const Real leaves_i_exp( const unsigned int i,
			     const Real alpha,
			     const Real r0 ) const;

    const Real p_leavea_i_exp( const unsigned int i,
			       const Real alpha,
			       const Real r0 ) const;

    const Real p_leaves_i_exp( const unsigned int i,
			       const Real alpha,
			       const Real r0 ) const;

    const Real p_int_r_i_exp( const unsigned int i,
			      const Real t,
			      const Real r,
			      const Real r0 ) const;

    const Real p_theta_n( const unsigned int n,
			  const RealVector& p_nTable, 
			  const RealVector& lgndTable ) const;

    const Real ip_theta_n( const unsigned int n,
			   const RealVector& p_nTable, 
			   const RealVector& lgndTable1 ) const;


    const Real p_int_r_i_exp_table( const unsigned int i,
				    const Real t,
				    const Real r,
				    const Real r0,
				    const RealVector& num_r0Table ) const;

    void initializeAlphaTable( const unsigned int n ) const;
    void updateAlphaTable0( const Real t ) const;
    void updateAlphaTable( const unsigned int n, const Real t ) const; 

    void createPsurvTable( RealVector& table, const Real r0 ) const; 
    void createNum_r0Table( RealVector& table, const Real r0 ) const;

    void createPleaveFactorTable( RealVector& table,
                                  const Real r0 ) const;
    void createPleavesTable( RealVector& table, const Real r0,
                             const RealVector& pleaveFactorTable ) const;
    void createPleaveaTable( RealVector& table, const Real r0,
                             const RealVector& pleaveFactorTable ) const;

    void makep_nTable( RealVector& p_nTable,
		       const Real r, 
		       const Real r0, 
		       const Real t ) const;
    
    void makedp_n_at_aTable( RealVector& p_nTable,
			     const Real r0, 
			     const Real t ) const;

    const Real findRoot( gsl_function& F,
                         gsl_root_fsolver* solver,
                         const Real low,
                         const Real high,
                         std::string funcName ) const;



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
	const Integer n;
	Real value;
    };

    static const Real 
    f_alpha_aux_F( const Real alpha,
		   const f_alpha_aux_params* const params );

    
    struct p_survival_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Real r0;
	const RealVector& psurvTable;
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
//	const RealVector& num_r0Table;
	const Real rnd;
    };

    static const Real 
    p_int_r_F( const Real r,
	       const p_int_r_params* const params );

    struct ip_theta_params
    { 
	const FirstPassagePairGreensFunction* const gf;
	const Real r;
	const Real r0;
	const Real t;
	const RealVector& p_nTable;
	const Real value;
    };

    static const Real 
    ip_theta_F( const Real theta,
		const ip_theta_params* const params );

    
    const Real num_r0( const Real alpha,
		       const Real r0 ) const;

    const Real pleaveFactor( const Real alpha,
                             const Real r0 ) const;

    static const Real P( const Integer n, const Real x );
    static const Real Q( const Integer n, const Real x );
    static const boost::tuple<Real,Real> P2( const Integer n, const Real x );
    static const boost::tuple<Real,Real> Q2( const Integer n, const Real x );



private:
    
    const Real h;
    const Real hsigma_p_1;

    mutable boost::array<Integer,MAX_ORDER+1> alphaOffsetTable;
    mutable boost::array<RealVector,MAX_ORDER+1> alphaTable;

    Real a;

};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION_HPP
