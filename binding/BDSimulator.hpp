#ifndef BINDING_BD_SIMULATOR_HPP
#define BINDING_BD_SIMULATOR_HPP

#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include "peer/util/to_native_converter.hpp"
#include "peer/wrappers/generator/generator_wrapper.hpp"
#include "peer/converters/tuple.hpp"
#include "peer/util/shared_const_ptr.hpp"

namespace binding {

template<typename Timpl>
void register_bd_simulator_class(char const* name)
{
    using namespace boost::python;
    using boost::shared_ptr;
    typedef Timpl impl_type;
    class_<impl_type, bases<typename impl_type::base_type>, boost::noncopyable>(
            name,
            init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&>())
        .def(init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&, double>())
        .def(init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&, double, int>())
//        .def("check", &impl_type::check)
//        .def("__len__", &impl_type::num_domains)
//        .def("__getitem__", &impl_type::get_domain)
//        .def("__iter__", &impl_type::get_domains,
//                return_value_policy<return_by_value>())
        ;

    peer::util::register_shared_const_ptr_from_python<typename impl_type::network_rules_type>();
}

} // namespace binding
#endif /* BINDING_BD_SIMULATOR_HPP */
