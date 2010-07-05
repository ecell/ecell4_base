#ifndef BINDING_PYEVENT_SCHEDULER_HPP
#define BINDING_PYEVENT_SCHEDULER_HPP

#include <boost/python.hpp>
#include "peer/converters/tuple.hpp"
#include "peer/converters/iterator.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base
register_event_scheduler_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    peer::converters::register_tuple_converter<
            typename impl_type::value_type>();

    peer::converters::register_stl_iterator_range_converter<
            typename impl_type::events_range, void*, return_by_value>();

    peer::converters::register_tuple_converter<
            typename impl_type::value_type>();

    return class_<impl_type, boost::noncopyable>(name)
        .add_property("time", &impl_type::time)
        .add_property("top",
            make_function(&impl_type::top,
                return_value_policy<return_by_value>()))
        .add_property("size", &impl_type::size)
        .add_property("second",
            make_function(&impl_type::second,
                return_value_policy<return_by_value>()))
        .def("update", &impl_type::update)
        .def("pop", &impl_type::pop,
            return_value_policy<return_by_value>())
        .def("clear", &impl_type::clear)
        .def("add", &impl_type::add)
        .def("check", &impl_type::check)
        .def("__getitem__", &impl_type::get,
            return_value_policy<return_by_value>())
        .def("__delitem__", &impl_type::remove)
        .def("__iter__", &impl_type::events,
            return_value_policy<return_by_value>())
        ;
}

} // namespace binding

#endif /* BINDING_PYEVENT_SCHEDULER_HPP */
