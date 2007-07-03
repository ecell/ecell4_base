#if !defined( __FREEGREENSFUNCTION )
#define __FREEGREENSFUNCTION 

#include <cmath>

#include <gsl/gsl_integration.h>

#include "GreensFunction.hpp"

/**
  Green's Function for a free diffusion particle.
*/

class FreeGreensFunction
    :
    public GreensFunction
{
    
private:

    static const Real TOLERANCE = 1e-8;

    static const Real H = 6;
    
public:
    
    FreeGreensFunction( const Real D )
        :
        GreensFunction( D )
    {
        ; // do nothing
    }
    
    
    ~FreeGreensFunction()
    {
        ; // do nothing
    }
    
    
    const Real drawTime( const Real ) const
    {
        return INFINITY;
    }
    
    const Real drawR( const Real rnd, 
		      const Real t ) const;
    
    const Real p_r( const Real r, 
                    const Real t ) const;

    const Real ip_r( const Real r, 
                     const Real t ) const;
    

private:

    struct ip_r_params
    { 
	const FreeGreensFunction* const gf;
	const Real t;
	const Real value;
    };


    static const Real ip_r_F( const Real r,
                              const ip_r_params* params );

};



#endif // __FREEGREENSFUNCTION 
