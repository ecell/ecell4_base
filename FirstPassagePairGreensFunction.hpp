#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION_HPP )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION_HPP 

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>

#include <gsl/gsl_roots.h>

#include "Logger.hpp"

#include "PairGreensFunction.hpp"


class FirstPassagePairGreensFunction: public PairGreensFunction
{
public:
    typedef std::vector<Real> RealVector;

private:
    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.
    static const Real THETA_TOLERANCE = 1e-5;

    static const Real MIN_T_FACTOR = 1e-8;

    static const unsigned int MAX_ORDER = 50;
    static const unsigned int MAX_ALPHA_SEQ = 2000;


public:
    
    FirstPassagePairGreensFunction(Real D, Real kf, Real Sigma, Real a);
    
    ~FirstPassagePairGreensFunction();

    Real geth() const
    {
        return this->h;
    }

    Real geta() const
    {
        return this->a;
    }
    
    Real drawTime(Real rnd, Real r0 ) const;

    std::pair<Real, EventType> 
    drawTime2(Real rnd1, Real rnd2, Real r0) const;

    EventType drawEventType(Real rnd, Real r0, Real t) const;
    
    Real drawR(Real rnd, Real r0, Real t) const;
    
    Real drawTheta(Real rnd, Real r, Real r0, Real t) const;
    
    Real f_alpha0(Real alpha) const;
    Real f_alpha0_aux(Real alpha) const;
  
    Real f_alpha(Real alpha, Integer n) const;
    Real f_alpha_aux(Real alpha, Integer n) const;

    Real p_0(Real t, Real r, Real r0) const;
    
    Real p_survival(Real t, Real r0) const;

    Real p_survival_table(Real t, Real r0, RealVector& table) const;

    Real p_leave_table(Real t, Real r0, RealVector const& table) const;


    Real dp_survival(Real t, Real r0) const;

    Real leaves(Real t, Real r0) const;

    Real leavea(Real t, Real r0) const;

    Real p_leaves(Real t, Real r0) const;

    Real p_leavea(Real t, Real r0) const;

    Real p_int_r(Real r, Real t, Real r0) const;

    Real p_theta(Real theta, Real r, Real r0, Real t) const;

    Real ip_theta(Real theta, Real r, Real r0, Real t) const;

    Real dp_theta(Real theta, Real r, Real r0, Real t) const;

    Real idp_theta(Real theta, Real r, Real r0, Real t) const;

    Real p_n(Integer n, Real r, Real r0, Real t, Real max_alpha) const;

    Real dp_n_at_a(Integer n, Real r0, Real t, Real max_alpha) const;


    Real p_n_alpha(unsigned int i, unsigned int n, Real r, Real r0, Real t) const;

    Real dp_n_alpha_at_a(unsigned int i, unsigned int n, Real r0, Real t) const;

    // methods below are kept public for debugging purpose.

    std::string dump() const;

    unsigned int alphaOffset(unsigned int n) const;

    Real alpha0_i(Integer i) const;

    Real alpha_i(Integer i, Integer n, gsl_root_fsolver* solver ) const;

    Real p_survival_i(Real alpha, Real r0) const;

    Real p_0_i(Real alpha, Real r, Real r0) const;

    Real dp_survival_i(Real alpha, Real r0) const;

    Real leavea_i(Real alpha, Real r0) const;

    Real leaves_i(Real alpha, Real r0) const;

    Real p_leavea_i(Real alpha, Real r0, Real pleave_factor) const;

    Real p_leaves_i(Real alpha, Real r0, Real pleave_factor) const;

    Real p_survival_den(Real alpha, Real r0) const;

    Real p_int_r_i(Real r, Real alpha, Real r0, Real num_r0) const;

    Real p_0_i_exp(unsigned int i, Real t, Real r, Real r0) const;

    Real p_survival_i_exp(unsigned int i, Real t, Real r0) const;

    Real p_survival_i_alpha(Real alpha, Real t, Real r0) const;


    Real p_survival_2i_exp(unsigned int i, Real t, Real r0) const;


protected:

    void clearAlphaTable() const;

    RealVector& getAlphaTable(size_t n) const
    {
        return this->alphaTable[n];
    }

    Real getAlpha(size_t n, RealVector::size_type i) const
    {
        RealVector& alphaTable( this->alphaTable[n] );
        RealVector::size_type oldSize( alphaTable.size() );

        if( oldSize <= i )
        {
            alphaTable.resize( i+1 );
            unsigned int offset( alphaOffset( n ) );

            gsl_root_fsolver* solver(
                gsl_root_fsolver_alloc(gsl_root_fsolver_brent));

            for( RealVector::size_type m( oldSize ); m <= i; ++m )
            {
                alphaTable[m] = alpha_i( m + offset, n, solver );
            }

            gsl_root_fsolver_free( solver );
        }

        return alphaTable[i];

    }

    Real getAlpha0(RealVector::size_type i) const
    {
        RealVector& alphaTable( this->alphaTable[0] );
        
        RealVector::size_type oldSize( alphaTable.size() );

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


    Real p_int_r_table(Real r, Real t, Real r0,
                       RealVector const& num_r0Table) const;

    Real ip_theta_table(Real theta, Real r, Real r0, Real t,
                        RealVector const& p_nTable) const;

    Real dp_theta_at_a(Real theta, Real r0, Real t ) const;


    Real p_theta_table(Real theta, Real r, Real r0, Real t, 
                       RealVector const& p_nTable ) const;

    void make_p_thetaTable( RealVector& pTable, Real r, Real r0, Real t,
                            unsigned int n, RealVector const& p_nTable ) const;

    Real p_survival_i_exp_table(unsigned int i, Real t, Real r0,
                                RealVector const& table ) const;

    Real p_leave_i_exp_table(unsigned int i, Real t, Real r0,
                             RealVector const& table ) const;


    Real dp_survival_i_exp(unsigned int i, Real alpha, Real r0) const;

    Real leavea_i_exp(unsigned int i, Real alpha, Real r0) const;

    Real leaves_i_exp(unsigned int i, Real alpha, Real r0) const;

    Real p_leavea_i_exp(unsigned int i, Real alpha, Real r0) const;

    Real p_leaves_i_exp(unsigned int i, Real alpha, Real r0) const;

    Real p_int_r_i_exp(unsigned int i, Real t, Real r, Real r0) const;

    Real p_int_r_i_exp_table(unsigned int i, Real t, Real r, Real r0,
                             RealVector& num_r0Table ) const;

    void initializeAlphaTable(unsigned int n) const;
    void updateAlphaTable0(Real t) const;
    void updateAlphaTable(unsigned int n, Real t) const; 

    void createPsurvTable(RealVector& table, Real r0) const; 
    void createNum_r0Table(RealVector& table, Real r0) const;

    void createPleaveFactorTable(RealVector& table, Real r0) const;
    void createPleavesTable(RealVector& table, Real r0,
                            RealVector const& pleaveFactorTable) const;
    void createPleaveaTable(RealVector& table, Real r0,
                            RealVector const& pleaveFactorTable) const;

    void makep_nTable(RealVector& p_nTable, Real r, Real r0, Real t) const;
    
    void makedp_n_at_aTable(RealVector& p_nTable, Real r0, Real t) const;

    unsigned int guess_maxi(Real t) const;

    Real 
    drawPleaves(gsl_function const& F,
                gsl_root_fsolver* solver,
                Real r0,
                Real t_guess,
                RealVector& pleaveFactorTable,
                RealVector& pleavesTable) const;

    Real 
    drawPleavea(gsl_function const& F,
                gsl_root_fsolver* solver,
                Real r0,
                Real t_guess,
                RealVector& pleaveFactorTable,
                RealVector& pleavesTable) const;

    
    Real num_r0(Real alpha, Real r0) const;

    Real pleaveFactor(Real alpha, Real r0) const;

    struct ip_theta_params;
    static Real ip_theta_F(Real, ip_theta_params const*);


private:
    
    const Real h;
    const Real hsigma_p_1;

    mutable boost::array<Integer, MAX_ORDER+1> alphaOffsetTable;
    mutable boost::array<RealVector, MAX_ORDER+1> alphaTable;

    const Real a;

    static Logger& log_;
};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION_HPP
