#ifndef BINDING_PYEVENT_HPP
#define BINDING_PYEVENT_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::object register_py_event_class(char const* name)
{
    using namespace boost::python;
    return class_<Timpl>(name, init<Real, object const&, object const&>())
        .def( "setTime", &PyEvent::setTime )
        .def( "getTime", &PyEvent::getTime )
        .def( "getObj", &PyEvent::getObj,
              return_value_policy<copy_const_reference>() )
        .def( "getArg", &PyEvent::getArg,
              return_value_policy<copy_const_reference>() );
}

} // namespace binding

#endif /* BINDING_PYEVENT_HPP */
