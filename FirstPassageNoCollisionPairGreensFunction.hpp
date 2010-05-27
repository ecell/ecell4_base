#if !defined( __FIRSTPASSAGENOCOLLISIONPAIRGREENSFUNCTION )
#define __FIRSTPASSAGENOCOLLISIONPAIRGREENSFUNCTION 

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>

#include <gsl/gsl_roots.h>

#include "Logger.hpp"
#include "PairGreensFunction.hpp"

class FirstPassageNoCollisionPairGreensFunction: public PairGreensFunction
{
public:
    typedef std::vector<Real> RealVector;

private:
    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.
    static const Real THETA_TOLERANCE = 1e-5;

    static const Real MIN_T = 1e-18;

    static const unsigned int MAX_ORDER = 50;
    static const unsigned int MAX_ALPHA_SEQ = 1005;


public:
    
    FirstPassageNoCollisionPairGreensFunction(Real D, Real a); 
    
    ~FirstPassageNoCollisionPairGreensFunction();

    Real geta() const
    {
        return this->a;
    }

    Real drawTime(Real rnd, Real r0) const;

    EventType drawEventType(Real rnd, Real r0, Real t) const;
    
    Real drawR(Real rnd, Real r0, Real t) const;
    
    Real drawTheta(Real rnd, Real r, Real r0, Real t) const;
    
    Real p_survival(Real t, Real r0) const;

    Real dp_survival(Real t, Real r0) const;

    Real p_int_r(Real r, Real t, Real r0) const;

    Real p_theta(Real theta, Real r, Real r0, Real t) const;

    Real ip_theta(Real theta, Real r, Real r0, Real t) const;

    Real dp_theta(Real theta, Real r, Real r0, Real t) const;

    Real idp_theta(Real theta, Real r, Real r0, Real t) const;


    Real p_n(Integer n, Real r, Real r0, Real t) const;

    Real dp_n(Integer n, Real r0, Real t ) const;


    Real p_n_alpha(unsigned int i, unsigned int n, Real r, Real r0, Real t ) const;

    Real dp_n_alpha(unsigned int i, unsigned int n, Real r0, Real t) const;

    // methods below are kept public for debugging purpose.

    std::string dump() const;

protected:

    Real p_theta_table(Real theta, Real r, Real r0, Real t, 
                       RealVector const& p_nTable ) const;

    Real ip_theta_table(Real theta, Real r, Real r0, Real t,
                        RealVector const& p_nTable ) const;

    void makep_nTable(RealVector& p_nTable,
                      Real r, Real r0, Real t) const;
    
    void makedp_nTable(RealVector& p_nTable, Real r0, Real t) const;

    struct ip_theta_params;
    static Real ip_theta_F(Real theta, ip_theta_params const* params);

private:
    
    mutable boost::array<Integer,MAX_ORDER+1> alphaOffsetTable;
    mutable boost::array<RealVector,MAX_ORDER+1> alphaTable;
    //mutable std::vector<RealVector> alphaTable;

    Real a;

    static Logger& log_;
};

#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION 
