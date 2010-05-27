#if !defined(__PLAINPAIRGREENSFUNCTION)
#define __PLAINPAIRGREENSFUNCTION 

#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>

#include "Logger.hpp"
#include "PairGreensFunction.hpp"

class BasicPairGreensFunction: public PairGreensFunction
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
    
    BasicPairGreensFunction(Real D, Real kf, Real Sigma);
    
    ~BasicPairGreensFunction();
    
    
    Real drawTime(Real rnd, Real r0) const;
    
    Real drawR(Real rnd, Real r0, Real t) const;
    
    Real drawTheta(Real rnd, Real r, Real r0, Real t) const;
    
    Real getkD() const
    {
        return this->kD;
    }
    
    Real getalpha() const
    {
        return this->alpha;
    }
    
    Real p_reaction(Real t, Real r0) const;
    Real p_survival(Real t, Real r0) const;
    Real p_int_r(Real r, Real t, Real r0) const;

//    const Real p_int_r_max(const Real t, const Real r0) const;

    
    Real p_theta(Real theta, Real r, Real r0, Real time) const;

    Real ip_theta(Real theta, Real r, Real r0, Real time) const;

    Real p_free(Real theta, Real r, Real r0, Real t) const;

    Real ip_free(Real theta, Real r, Real r0, Real t) const;
    
    Real p_corr(Real theta, Real r, Real r0, Real t) const;

    Real ip_corr(Real theta, Real r, Real r0, Real t) const;

    std::string dump() const;

private:
    Real p_corr_R(Real alpha, unsigned int n, Real r, Real r0, Real t) const;

    
    Real p_corr_n(unsigned int n, RealVector const& RnTable, RealVector const& lgndTable) const;

    Real ip_corr_n(unsigned int n, RealVector const& RnTable, RealVector const& lgndTable) const;

    Real p_corr_table(Real theta, Real r, Real r0, Real t, RealVector const& RnTable) const;

    Real 
    ip_corr_table(Real theta, Real r, Real r0, Real t, RealVector const& RnTable) const;
    
    Real p_theta_table(Real r, Real r0, Real theta, Real time,
                       RealVector const& RnTable) const;

    Real ip_theta_table(Real r, Real r0, Real theta, Real time,
                        RealVector const& RnTable) const;

    Real 
    p_corr_table(Real theta, Real r, Real r0, Real t, RealVector const& RnTable);
    

    void makeRnTable(RealVector& RnTable, Real r, Real r0, Real t) const;

    Real Rn(unsigned int order, Real r, Real r0, Real t,
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
