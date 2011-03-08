#if !defined( __FIRSTPASSAGENOCOLLISIONPAIRGREENSFUNCTION )
#define __FIRSTPASSAGENOCOLLISIONPAIRGREENSFUNCTION 

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>

#include <gsl/gsl_roots.h>

#include "Logger.hpp"
#include "GreensFunction3DRadAbsBase.hpp"

class GreensFunction3DAbs: public GreensFunction3DRadAbsBase
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
    
    GreensFunction3DAbs(Real D, Real r0, Real a); 
    
    virtual ~GreensFunction3DAbs();

    Real geta() const
    {
        return this->a;
    }

    virtual Real drawTime(Real rnd) const;

    virtual EventKind drawEventType(Real rnd, Real t) const;
    
    virtual Real drawR(Real rnd, Real t) const;
    
    virtual Real drawTheta(Real rnd, Real r, Real t) const;
    
    Real p_survival(Real t) const;

    Real dp_survival(Real t) const;

    Real p_int_r(Real r, Real t) const;

    Real p_theta(Real theta, Real r, Real t) const;

    Real ip_theta(Real theta, Real r, Real t) const;

    Real dp_theta(Real theta, Real r, Real t) const;

    Real idp_theta(Real theta, Real r, Real t) const;


    Real p_n(Integer n, Real r, Real t) const;

    Real dp_n(Integer n, Real t ) const;


    Real p_n_alpha(unsigned int i, unsigned int n, Real r, Real t ) const;

    Real dp_n_alpha(unsigned int i, unsigned int n, Real t) const;

    // methods below are kept public for debugging purpose.

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction3DAbs";
    }

protected:

    Real p_theta_table(Real theta, Real r, Real t, 
                       RealVector const& p_nTable ) const;

    Real ip_theta_table(Real theta, Real r, Real t,
                        RealVector const& p_nTable ) const;

    void makep_nTable(RealVector& p_nTable, Real r, Real t) const;
    
    void makedp_nTable(RealVector& p_nTable, Real t) const;

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
