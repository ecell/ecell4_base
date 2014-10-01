#if !defined( __FREEGREENSFUNCTION )
#define __FREEGREENSFUNCTION 

#include "compat.h"

#include <gsl/gsl_integration.h>

#include "Logger.hpp"
#include "GreensFunction.hpp"

/**
  Green's Function for a free diffusion particle.
*/

class GreensFunction3DSym
    :
    public GreensFunction
{
    
private:

    static const Real TOLERANCE = 1e-8;

    static const Real H = 6;
    
public:
    
    GreensFunction3DSym( const Real D )
        :
        GreensFunction( D )
    {
        ; // do nothing
    }
    
    
    ~GreensFunction3DSym()
    {
        ; // do nothing
    }

    Real drawTime( const Real ) const
    {
        return INFINITY;
    }
    
    Real drawR( const Real rnd, const Real t ) const;
    
    Real p_r( const Real r, const Real t ) const;

    Real ip_r( const Real r, const Real t ) const;
    

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction3DSym";
    }

private:

    static Logger& log_;
};



#endif // __FREEGREENSFUNCTION 
