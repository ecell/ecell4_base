#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION 

#include <iostream>
#include <math.h>

#include <gsl/gsl_integration.h>

#include "PairGreensFunction.hpp"

class FirstPassagePairGreensFunction
    :
    public PairGreensFunction
{
    
    
public:
    
    FirstPassagePairGreensFunction( const Real D, 
				    const Real kf, 
				    const Real Sigma );
    
    virtual ~FirstPassagePairGreensFunction();
    
    
    virtual const Real drawTime( const Real rnd, const Real r0,
				 const Real maxt ) const;
    
    virtual const Real drawR( const Real rnd, 
			      const Real r0, 
			      const Real t ) const;
    
    virtual const Real drawTheta( const Real rnd,
				  const Real r, 
				  const Real r0, 
				  const Real t ) const;
    
    
    const Real getkD() const
    {
	return this->kD;
    }
    
    const Real getalpha() const
    {
	return this->alpha;
    }
    

  

private:
    
    
    
private:
    
    const Real kD;
    const Real alpha;
    
    static const Real P_CUTOFF = 1e-6;
    
};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION 
