#ifndef BINDING_PYEVENT_HPP
#define BINDING_PYEVENT_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::object register_event_class(char const* name)
{
    using namespace boost::python;
    return class_<Timpl>(name, init<Real, object const&>())
        .add_property("time",
            make_function(
                &Timpl::time, return_value_policy<copy_const_reference>()))
        .add_property("data",
            make_function(&Timpl::data,
                          return_value_policy<copy_const_reference>()))
        ;
}

} // namespace binding

#endif /* BINDING_PYEVENT_HPP */
