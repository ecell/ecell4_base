#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <gsl/gsl_errno.h>

#include "findRoot.hpp"



const Real 
findRoot( gsl_function& F,
          gsl_root_fsolver* solver,
          const Real low,
          const Real high,
          const Real tol_abs,
          const Real tol_rel,
          std::string funcName )
{
    Real l( low );
    Real h( high );

    gsl_root_fsolver_set( solver, &F, l, h );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
        gsl_root_fsolver_iterate( solver );
        l = gsl_root_fsolver_x_lower( solver );
        h = gsl_root_fsolver_x_upper( solver );

        const int status( gsl_root_test_interval( l, h, tol_abs,
                                                  tol_rel ) );

        if( status == GSL_CONTINUE )
        {
            if( i >= maxIter )
            {
                gsl_root_fsolver_free( solver );
                std::cerr << funcName << ": failed to converge." << std::endl;
                throw std::exception();
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  

    const Real root( gsl_root_fsolver_root( solver ) );

    return root;
}
