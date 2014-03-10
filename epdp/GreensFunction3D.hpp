#if !defined( __FREEPAIRGREENSFUNCTION )
#define __FREEPAIRGREENSFUNCTION 

#include "compat.h"

#include <gsl/gsl_integration.h>

#include "Logger.hpp"
#include "PairGreensFunction.hpp"

/**
   Pair Green's function for the case where the pair never interact.

   Therefore, drawTime() always returns +INFINITY.
   kf == sigma == 0.
*/

class GreensFunction3D: public PairGreensFunction
{
    
private:

    static const Real TOLERANCE = 1e-8;

    static const Real H = 7;
    
public:
    
    GreensFunction3D(Real D, Real r0)
        : PairGreensFunction(D, 0.0, r0, 0.0)
    {
        ; // do nothing
    }
    
    
    virtual ~GreensFunction3D();
    
    virtual Real drawTime(Real rnd) const;

    Real drawR(Real rnd, Real t) const;
    
    Real drawTheta(Real rnd, Real r, Real t) const;

    Real p_r(Real r, Real t) const;

    Real ip_r(Real r, Real t ) const;
    

    Real p_theta(Real theta, Real r, Real t) const;

    Real ip_theta(Real theta, Real r, Real t ) const;

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction3D";
    }

private:
    static Logger& log_;
};


#endif // __PLAINPAIRGREENSFUNCTION 
