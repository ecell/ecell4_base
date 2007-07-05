#include <boost/python.hpp>

#include "EventScheduler.hpp"

//using namespace boost::python;


class PyEvent
    :
    public libecs::EventBase
{
    
public:
	
    PyEvent( const double time, const boost::python::object& obj )
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
    

    const boost::python::object& getObj() const
    {
	return this->obj;
    }

    void fire()
    {
	boost::python::object ret( this->obj.attr( "fire" )() );
	this->setTime( this->getTime() + 
                       boost::python::extract<double>( ret ) );
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

    boost::python::object obj;
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
    
    const EventID addEvent( const double t, 
			    const boost::python::object& obj )
    {
	return libecs::EventScheduler<PyEvent>::addEvent( PyEvent( t, obj ) );
    }


    void updateEvent( const EventID id, const double t, 
                      const boost::python::object& obj )
    {
	libecs::EventScheduler<PyEvent>::updateEvent( id,
                                                      PyEvent( t, obj ) );
    }


};



/*
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
*/
