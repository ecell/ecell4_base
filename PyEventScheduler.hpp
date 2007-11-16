#include <boost/python.hpp>

#include "EventScheduler.hpp"

//using namespace boost::python;


class PyEvent
    :
    public libecs::EventBase
{
    
public:
	
    PyEvent( const double time, 
             const boost::python::object& obj,
             const boost::python::object& arg )
        :
        EventBase( time ),
        obj( obj ),
        arg( arg )
    {
	; // do nothing
    }

    ~PyEvent()
    {
	; // do nothing
    }
    

    const boost::python::object& getObj() const
    {
	return this->obj;
    }

    const boost::python::object& getArg() const
    {
	return this->arg;
    }

    void fire()
    {
	boost::python::object ret( this->obj( this->arg ) );
	this->setTime( this->getTime() + 
                       boost::python::extract<double>( ret ) );
    }

    PyEvent() // dummy
    {
	; // do nothing
    }

private:

    boost::python::object obj;
    boost::python::object arg;
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
			    const boost::python::object& func,
			    const boost::python::object& arg )
    {
	return libecs::EventScheduler<PyEvent>::
            addEvent( PyEvent( t, func, arg ) );
    }


    void updateEventTime( const EventID id, const double t ) 
    {
	libecs::EventScheduler<PyEvent>::updateEventTime( id, t );
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
