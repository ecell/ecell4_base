#include <exception>

#include <gsl/gsl_math.h>

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <numpy/arrayobject.h>

#include "DynamicPriorityQueue.hpp"
#include "EventScheduler.hpp"

using namespace boost::python;


class PyEvent
    :
    public libecs::EventBase
{
    
public:
	
    PyEvent( const double time, object obj )
        :
        EventBase( time ),
	obj( obj )
    {
	; // do nothing
    }

    virtual ~PyEvent()
    {
	; // do nothing
    }
    

    const object& getObj() const
    {
	return this->obj;
    }

    void fire()
    {
	object ret( this->obj.attr( "fire" )() );
	this->setTime( this->getTime() + extract<double>(ret) );
    }

    void update( const double t )
    {
	this->obj.attr( "update" )( t );
    }

    const bool isDependentOn( const PyEvent& arg ) const
    {
	return this->obj.attr( "isDependentOn" )( arg.getObj() );
    }

    PyEvent() // dummy
    {
	; // do nothing
    }

private:

    object obj;
};



class PyEventScheduler
    :
    public libecs::EventScheduler<PyEvent>
{
public:

    typedef libecs::EventScheduler<PyEvent>::EventIndex EventIndex;
    
    PyEventScheduler()
    {
	; // do nothing
    }
    
    ~PyEventScheduler()
    {
	; // do nothing
    }
    
    const EventIndex addEvent( const double t, const object& obj )
    {
	return libecs::EventScheduler<PyEvent>::addEvent( PyEvent( t, obj ) );
    }


};




class PyEvent_to_python
{
public:

  static PyObject* 
  convert( const PyEvent& value )
  {
      return PyTuple_Pack( 2, 
			   PyFloat_FromDouble( value.getTime() ),
			   value.getObj().ptr() );
  }

};






#include "PairGreensFunction.hpp"
#include "PlainPairGreensFunction.hpp"
#include "FirstPassageGreensFunction.hpp"
#include "FirstPassagePairGreensFunction.hpp"

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


BOOST_PYTHON_MODULE( _gfrd )
{
    import_array();

    to_python_converter< PyEvent, PyEvent_to_python>();
  

    class_<PyEvent, boost::noncopyable>( "Event", init<const Real,
					 object>() )
	.def( "setTime", &PyEvent::setTime )
	.def( "getTime", &PyEvent::getTime )
//	.def( "fire", &PyEvent::fire )
//	.def( "update", &PyEvent::update )
//	.def( "isDependentOn", &PyEvent::isDependentOn )
	;

    typedef const PyEventScheduler::Event& 
	(PyEventScheduler::*geteventrefsig)() const;
    typedef PyEventScheduler::Event& (PyEventScheduler::*geteventrefbyindexsig)
	( const PyEventScheduler::EventIndex );

    class_<PyEventScheduler, boost::noncopyable>( "EventScheduler" )
	.def( "getTime", &PyEventScheduler::getTime )
	.def( "getSize", &PyEventScheduler::getSize )
	.def( "getTopEvent", geteventrefsig( &PyEventScheduler::getTopEvent ),
//	      return_value_policy<reference_existing_object>() )
//	      return_value_policy<copy_const_reference>() )
	      return_value_policy<return_by_value>() )
	.def( "getTopIndex", &PyEventScheduler::getTopIndex )
	.def( "getEvent", geteventrefbyindexsig( &PyEventScheduler::getEvent ),
	      return_value_policy<reference_existing_object>() )
	.def( "step", &PyEventScheduler::step )
	.def( "clear", &PyEventScheduler::clear )
	.def( "addEvent", &PyEventScheduler::addEvent )
	.def( "updateAllEventDependency", 
	      &PyEventScheduler::updateAllEventDependency )
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

    class_<FirstPassageGreensFunction>( "FirstPassageGreensFunction",
					init<const Real>() )
	.def( "getD", &FirstPassageGreensFunction::getD )
	.def( "drawTime", &FirstPassageGreensFunction::drawTime )
	.def( "drawR", &FirstPassageGreensFunction::drawR )
	.def( "p_survival", &FirstPassageGreensFunction::p_survival )
	.def( "p_r_int", &FirstPassageGreensFunction::p_r_int )
	.def( "p_r_fourier", &FirstPassageGreensFunction::p_r_fourier )
	;


    class_<FirstPassagePairGreensFunction>( "FirstPassagePairGreensFunction",
					    init<const Real, 
					    const Real,
					    const Real>() )
	.def( "getD", &FirstPassagePairGreensFunction::getD )
	.def( "drawTime", &FirstPassagePairGreensFunction::drawTime )
	.def( "drawR", &FirstPassagePairGreensFunction::drawR )

	.def( "p_survival", &FirstPassagePairGreensFunction::p_survival )

	// debug
	.def( "f_alpha", &FirstPassagePairGreensFunction::f_alpha )
	.def( "f_alpha_survival", 
	      &FirstPassagePairGreensFunction::f_alpha_survival )
	.def( "alpha_survival_n", 
	      &FirstPassagePairGreensFunction::alpha_survival_n )
	;

    def( "distanceSq", &distanceSq_ );

}
