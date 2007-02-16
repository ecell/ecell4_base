#define HAVE_INLINE

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <exception>
#include <vector>

#include <boost/array.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
//#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>

#include "bessel.hpp"

#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassagePairGreensFunction.hpp"



FirstPassagePairGreensFunction::
FirstPassagePairGreensFunction( const Real D, 
				const Real kf, 
				const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
    kD( 4.0 * M_PI * getSigma() * getD() ),
    alpha( ( 1.0 + ( getkf() / getkD() ) ) * ( sqrt( getD() ) / getSigma() ) )
{
    ; // do nothing
}

FirstPassagePairGreensFunction::~FirstPassagePairGreensFunction()
{
    ; // do nothing
}




