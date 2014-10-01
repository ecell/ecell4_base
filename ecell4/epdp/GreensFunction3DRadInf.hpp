#if !defined(__PLAINPAIRGREENSFUNCTION)
#define __PLAINPAIRGREENSFUNCTION 

#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>

#include "Logger.hpp"
#include "PairGreensFunction.hpp"

class GreensFunction3DRadInf: public PairGreensFunction
{
public:
    typedef std::vector<Real> RealVector;

private:
    struct p_corr_R_params;
    struct p_theta_params;

private:
    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.
    static const Real THETA_TOLERANCE = 1e-5;

    static const Real MIN_T = 1e-12;

    static const unsigned int MAX_ORDER = 70;

    static const Real H = 4.0;
    

    
public:
    
    GreensFunction3DRadInf(Real D, Real kf, Real r0, Real Sigma);

    virtual ~GreensFunction3DRadInf();
 
    virtual Real drawTime(Real rnd) const;
    
    virtual Real drawR(Real rnd, Real t) const;
    
    virtual Real drawTheta(Real rnd, Real r, Real t) const;
    
    Real getkD() const
    {
        return this->kD;
    }
    
    Real getalpha() const
    {
        return this->alpha;
    }
    
    Real p_reaction(Real t) const;
    Real p_survival(Real t) const;
    Real p_int_r(Real r, Real t) const;
    
    Real p_theta(Real theta, Real r, Real time) const;

    Real ip_theta(Real theta, Real r, Real time) const;

    Real p_free(Real theta, Real r, Real t) const;

    Real ip_free(Real theta, Real r, Real t) const;
    
    Real p_corr(Real theta, Real r, Real t) const;

    Real ip_corr(Real theta, Real r, Real t) const;

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction3DRadInf";
    }

private:
    Real p_corr_R(Real alpha, unsigned int n, Real r, Real t) const;

    
    Real p_corr_n(unsigned int n, RealVector const& RnTable, RealVector const& lgndTable) const;

    Real ip_corr_n(unsigned int n, RealVector const& RnTable, RealVector const& lgndTable) const;

    Real p_corr_table(Real theta, Real r, Real t, RealVector const& RnTable) const;

    Real 
    ip_corr_table(Real theta, Real r, Real t, RealVector const& RnTable) const;
    
    Real p_theta_table(Real r, Real theta, Real time,
                       RealVector const& RnTable) const;

    Real ip_theta_table(Real r, Real theta, Real time,
                        RealVector const& RnTable) const;

    Real 
    p_corr_table(Real theta, Real r, Real t, RealVector const& RnTable);
    

    void makeRnTable(RealVector& RnTable, Real r, Real t) const;

    Real Rn(unsigned int order, Real r, Real t,
            gsl_integration_workspace* workspace, Real tol) const;

private:
    static Real p_corr_R_F(Real, p_corr_R_params*);
    static Real ip_theta_F(Real theta, p_theta_params* params);
    
private:
    const Real kD;
    const Real alpha;
   
    static Logger& log_;
};



#endif // __PLAINPAIRGREENSFUNCTION 
