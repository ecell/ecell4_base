//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <boost/bind.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sum.h>

#include "factorial.hpp"
#include "bessel.hpp"

//#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassageNoCollisionPairGreensFunction.hpp"



FirstPassageNoCollisionPairGreensFunction::
FirstPassageNoCollisionPairGreensFunction( const Real D ) 
    :
    PairGreensFunction( D, 0, 0 ),
    a( INFINITY )
{
    ; // do nothing
}

FirstPassageNoCollisionPairGreensFunction::
~FirstPassageNoCollisionPairGreensFunction()
{
    ; // do nothing
}

void FirstPassageNoCollisionPairGreensFunction::seta( const Real a )
{
    this->a = a;

//    clearAlphaTable();
}
