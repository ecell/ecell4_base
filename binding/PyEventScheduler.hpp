#ifndef BINDING_PYEVENT_SCHEDULER_HPP
#define BINDING_PYEVENT_SCHEDULER_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base
register_py_event_scheduler_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;
    typedef const typename impl_type::Event& (impl_type::*geteventrefsig)() const;
    typedef const typename impl_type::Event& (impl_type::*geteventrefbyindexsig)(
            const typename impl_type::EventID);

    return class_<impl_type, boost::noncopyable>(name)
        .def("getTime", &impl_type::getTime)
        .def("getTopTime", &impl_type::getTopTime)
        .def("getSize", &impl_type::getSize)
        .def("getTopEvent", geteventrefsig( &impl_type::getTopEvent),
              return_value_policy<copy_const_reference>())
        .def("getTopID", &impl_type::getTopID)
        .def("peekSecondEvent", 
              geteventrefsig( &impl_type::peekSecondEvent),
              return_value_policy<copy_const_reference>())
        .def("getEvent", geteventrefbyindexsig( &impl_type::getEvent),
              return_value_policy<copy_const_reference>())
        .def("getEventByIndex", &impl_type::getEventByIndex,
              return_value_policy<copy_const_reference>())
        .def("step", &impl_type::step)
        .def("clear", &impl_type::clear)
        .def("addEvent", &impl_type::addEvent)
        .def("removeEvent", &impl_type::removeEvent)
        .def("updateEventTime", &impl_type::updateEventTime)
//        .def("updateAllEventDependency", 
//              &impl_type::updateAllEventDependency)
        .def("check", &impl_type::check)
        ;
}

} // namespace binding

#endif /* BINDING_PYEVENT_SCHEDULER_HPP */
