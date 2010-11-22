#ifndef BINDING_PYTHON_EVENT_HPP
#define BINDING_PYTHON_EVENT_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::object register_python_event_class(char const* name)
{
    using namespace boost::python;

    return class_<Timpl, bases<typename Timpl::base_type>, boost::shared_ptr<Timpl>, boost::noncopyable>(name, init<typename Timpl::time_type>())
        .def(init<typename Timpl::time_type, boost::python::object const&>())
        .add_property("data",
            make_function(&Timpl::data,
                          return_value_policy<copy_const_reference>()))
        ;
}

} // namespace binding

#endif /* BINDING_PYEVENT_HPP */
