#if !defined( __FREEPAIRGREENSFUNCTION )
#define __FREEPAIRGREENSFUNCTION 

#include <cmath>

#include "compat.h"

#include <gsl/gsl_integration.h>

#include "Logger.hpp"
#include "PairGreensFunction.hpp"

/**
   Pair Green's function for the case where the pair never interact.

   Therefore, drawTime() always returns +INFINITY.
   kf == sigma == 0.
*/

class FreePairGreensFunction
    :
    public PairGreensFunction
{
    
private:

    static const Real TOLERANCE = 1e-8;

    static const Real H = 7;
    
public:
    
    FreePairGreensFunction(Real D, Real r0)
        : PairGreensFunction(D, 0.0, r0, 0.0)
    {
        ; // do nothing
    }
    
    
    ~FreePairGreensFunction()
    {
        ; // do nothing
    }
    
    
    Real drawTime(Real rnd, Real maxt) const
    {
        return INFINITY;
    }
    
    Real drawR(Real rnd, Real t) const;
    
    Real drawTheta(Real rnd, Real r, Real t) const;

    Real p_r(Real r, Real t) const;

    Real ip_r(Real r, Real t ) const;
    

    Real p_theta(Real theta, Real r, Real t) const;

    Real ip_theta(Real theta, Real r, Real t ) const;

    std::string dump() const;

private:
    static Logger& log_;
};


#endif // __PLAINPAIRGREENSFUNCTION 
