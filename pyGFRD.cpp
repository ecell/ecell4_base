#include <exception>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <boost/lexical_cast.hpp>
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <numpy/arrayobject.h>

#include "PyEventScheduler.hpp"

#include "PlainPairGreensFunction.hpp"
#include "FreePairGreensFunction.hpp"
#include "FirstPassageGreensFunction.hpp"
#include "FirstPassagePairGreensFunction.hpp"


using namespace boost::python;


const double distanceSq( const double* const p1, const double* const p2 )
{
    return gsl_pow_2( p1[0] - p2[0] ) 
	+ gsl_pow_2( p1[1] - p2[1] ) 
	+ gsl_pow_2( p1[2] - p2[2] );
}

const double distanceSq_( PyArrayObject* o1, PyArrayObject* o2 )
{
    // check type, dimension, size.
    if( o1->descr->type_num != PyArray_DOUBLE || 
	o2->descr->type_num != PyArray_DOUBLE ||
	o1->nd != 1 || o2->nd != 1 ||
	o1->dimensions[0] != 3 || o2->dimensions[0] != 3 )
    {
	throw std::exception();
    }

    const double* const 
	p1( reinterpret_cast<const double* const>( o1->data ) );
    const double* const 
	p2( reinterpret_cast<const double* const>( o2->data ) );
  
    return distanceSq( p1, p2 );
}


void* extract_pyarray(PyObject* x)
{
    return PyObject_TypeCheck(x, &PyArray_Type) ? x : 0;
}



// Exception translators here.

void translateException( const std::exception& anException )
{
  PyErr_SetString( PyExc_RuntimeError, anException.what() );
}


// GSL error handler.
void gfrd_gsl_error_handler( const char* reason,
			     const char* file,
			     int line,
			     int gsl_errno )
{
    throw std::runtime_error( std::string( "GSL error: " ) +
			      std::string( reason ) +
			      std::string( " at " ) +
			      std::string( file ) + std::string( ":" ) +
			      boost::lexical_cast<std::string>( line ) );
}


BOOST_PYTHON_MODULE( _gfrd )
{
    import_array();

    // GSL error handler: is this the best place for this?
    gsl_set_error_handler( &gfrd_gsl_error_handler );

  
    register_exception_translator<std::exception>( &translateException );


    to_python_converter<PyEvent, PyEvent_to_python>();

    class_<PyEvent, boost::noncopyable>( "Event", init<const Real,
					 const object&>() )
	.def( "setTime", &PyEvent::setTime )
	.def( "getTime", &PyEvent::getTime )
//	.def( "getObj", &PyEvent::getObj )
//	.def( "fire", &PyEvent::fire )
//	.def( "update", &PyEvent::update )
//	.def( "isDependentOn", &PyEvent::isDependentOn )
	;


    typedef const PyEventScheduler::Event& 
	(PyEventScheduler::*geteventrefsig)() const;
    typedef PyEventScheduler::Event& (PyEventScheduler::*geteventrefbyindexsig)
	( const PyEventScheduler::EventID );

    class_<PyEventScheduler, boost::noncopyable>( "EventScheduler" )
	.def( "getTime", &PyEventScheduler::getTime )
	.def( "getNextTime", &PyEventScheduler::getNextTime )
	.def( "getSize", &PyEventScheduler::getSize )
	.def( "getTopEvent", geteventrefsig( &PyEventScheduler::getTopEvent ),
//	      return_value_policy<reference_existing_object>() )
//	      return_value_policy<copy_const_reference>() )
	      return_value_policy<return_by_value>() )
	.def( "getEvent", geteventrefbyindexsig( &PyEventScheduler::getEvent ),
	      return_value_policy<return_by_value>() )
	.def( "getEventByIndex", &PyEventScheduler::getEventByIndex,
	      return_value_policy<return_by_value>() )
	.def( "step", &PyEventScheduler::step )
	.def( "clear", &PyEventScheduler::clear )
	.def( "addEvent", &PyEventScheduler::addEvent )
	.def( "removeEvent", &PyEventScheduler::removeEvent )
	.def( "updateEvent", &PyEventScheduler::updateEvent )
//	.def( "updateAllEventDependency", 
//	      &PyEventScheduler::updateAllEventDependency )
	;


    //    converter::registry::insert( &extract_pyarray, type_id<PyArrayObject>());


    class_<PlainPairGreensFunction>( "PlainPairGreensFunction",
				     init<const Real, 
				     const Real, 
				     const Real>() )
	.def( "getD", &PlainPairGreensFunction::getD )
	.def( "getkf", &PlainPairGreensFunction::getkf )
	.def( "getSigma", &PlainPairGreensFunction::getSigma )
	.def( "drawTime", &PlainPairGreensFunction::drawTime )
	.def( "drawR", &PlainPairGreensFunction::drawR )
	.def( "drawTheta", &PlainPairGreensFunction::drawTheta )
	;

    class_<FreePairGreensFunction>( "FreePairGreensFunction",
                                    init<const Real>() )
	.def( "getD", &FreePairGreensFunction::getD )
	.def( "getkf", &FreePairGreensFunction::getkf )
	.def( "getSigma", &FreePairGreensFunction::getSigma )
	.def( "drawTime", &FreePairGreensFunction::drawTime )
	.def( "drawR", &FreePairGreensFunction::drawR )
	.def( "drawTheta", &FreePairGreensFunction::drawTheta )

	.def( "p_theta", &FreePairGreensFunction::p_theta )
	.def( "ip_theta", &FreePairGreensFunction::ip_theta )
	;

    class_<FirstPassageGreensFunction>( "FirstPassageGreensFunction",
					init<const Real>() )
	.def( "getD", &FirstPassageGreensFunction::getD )
	.def( "drawTime", &FirstPassageGreensFunction::drawTime )
	.def( "drawR", &FirstPassageGreensFunction::drawR )
	.def( "p_survival", &FirstPassageGreensFunction::p_survival )
	.def( "p_r_int", &FirstPassageGreensFunction::p_r_int )
	.def( "p_r_fourier", &FirstPassageGreensFunction::p_r_fourier )
	;

    enum_<EventType>( "EventType" )
	.value( "REACTION", REACTION )
	.value( "ESCAPE", ESCAPE )
	;

    class_<FirstPassagePairGreensFunction>( "FirstPassagePairGreensFunction",
					    init<const Real, 
					    const Real,
					    const Real>() )
	.def( "seta", &FirstPassagePairGreensFunction::seta )
	.def( "geta", &FirstPassagePairGreensFunction::geta )
	.def( "getD", &FirstPassagePairGreensFunction::getD )
	.def( "getkf", &PlainPairGreensFunction::getkf )
	.def( "getSigma", &PlainPairGreensFunction::getSigma )
	.def( "drawTime", &FirstPassagePairGreensFunction::drawTime )
	.def( "drawEventType", &FirstPassagePairGreensFunction::drawEventType )
	.def( "drawR", &FirstPassagePairGreensFunction::drawR )
	.def( "drawTheta", &FirstPassagePairGreensFunction::drawTheta )

	.def( "p_survival", &FirstPassagePairGreensFunction::p_survival )
	.def( "dp_survival", &FirstPassagePairGreensFunction::dp_survival )
	.def( "p_leaves", &FirstPassagePairGreensFunction::p_leaves )
	.def( "p_leavea", &FirstPassagePairGreensFunction::p_leavea )
	.def( "leaves", &FirstPassagePairGreensFunction::leaves )
	.def( "leavea", &FirstPassagePairGreensFunction::leavea )
	.def( "p_0", &FirstPassagePairGreensFunction::p_0 )
	.def( "p_int_r", &FirstPassagePairGreensFunction::p_int_r )
	.def( "p_int_r", &FirstPassagePairGreensFunction::p_int_r )
	.def( "p_theta", &FirstPassagePairGreensFunction::p_theta )
	.def( "ip_theta", &FirstPassagePairGreensFunction::ip_theta )
	.def( "idp_theta", &FirstPassagePairGreensFunction::idp_theta )

	.def( "f_alpha0", &FirstPassagePairGreensFunction::f_alpha0 )
	.def( "alpha0_i", &FirstPassagePairGreensFunction::alpha0_i )
	.def( "f_alpha", &FirstPassagePairGreensFunction::f_alpha )
	.def( "f_alpha_aux", &FirstPassagePairGreensFunction::f_alpha_aux )

	.def( "dump", &FirstPassagePairGreensFunction::dump )

//	.def( "alpha_i", &FirstPassagePairGreensFunction::alpha_i )
	;

    def( "distanceSq", &distanceSq_ );

}
