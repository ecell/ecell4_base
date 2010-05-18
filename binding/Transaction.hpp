#ifndef BINDING_TRANSACTION_HPP
#define BINDING_TRANSACTION_HPP

#include <boost/python.hpp>
#include "utils/range.hpp"

namespace binding {

template<typename Timpl_, typename Tbase_>
inline boost::python::objects::class_base
register_transaction_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    return class_<impl_type, bases<Tbase_>, boost::noncopyable>(name, no_init)
        .add_property("added_particles",
            make_function(&impl_type::get_added_particles,
                return_value_policy<return_by_value>()))
        .add_property("removed_particles",
            make_function(&impl_type::get_removed_particles,
                return_value_policy<return_by_value>()))
        .add_property("modified_particles",
            make_function(&impl_type::get_modified_particles,
                return_value_policy<return_by_value>()))
        .def("rollback", &impl_type::rollback)
        ;
}

template<typename Timpl_, typename Tbase_>
inline boost::python::objects::class_base
register_transaction_impl_class(char const* name)
{
    using namespace boost::python;

    return class_<Timpl_, bases<Tbase_>, boost::noncopyable>(
            name, init<typename Timpl_::particle_container_type&>());
}

} // namesapce binding

#endif /* BINDING_TRANSACTION_HPP */
