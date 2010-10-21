#ifndef EGFRD_SIMULATOR_HPP
#define EGFRD_SIMULATOR_HPP

#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include "peer/util/to_native_converter.hpp"
#include "peer/wrappers/generator/generator_wrapper.hpp"
#include "peer/converters/tuple.hpp"

namespace binding {

template<typename Timpl_>
struct shell_variant_converter
{
    typedef Timpl_ impl_type;
    typedef typename impl_type::shell_variant_type native_type;
    typedef typename impl_type::spherical_shell_type spherical_shell_type;
    typedef typename impl_type::cylindrical_shell_type cylindrical_shell_type;

    struct visitor: boost::static_visitor<boost::python::object>
    {
        boost::python::object operator()(boost::none_t const&) const
        {
            return boost::python::object();
        }

        boost::python::object operator()(spherical_shell_type const& shell) const
        {
            return boost::python::object(shell);
        }

        boost::python::object operator()(cylindrical_shell_type const& shell) const
        {
            return boost::python::object(shell);
        }
    };

    static PyObject* convert(native_type const& val)
    {
        return boost::python::incref(boost::apply_visitor(visitor(), val).ptr());
    }

    static void __register()
    {
        boost::python::to_python_converter<native_type, shell_variant_converter>();
    }
};


template<typename Timpl>
void register_egfrd_simulator_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;
    typedef std::pair<typename impl_type::shell_id_type, typename impl_type::shell_variant_type> get_shell_result_type;
    class_<impl_type, bases<typename impl_type::base_type>, boost::noncopyable>(
            name,
            init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&>())
        .def(init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&, int>())
        .def(init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&, int,
                 typename impl_type::length_type>())
        .def("get_shell", static_cast<get_shell_result_type(impl_type::*)(typename impl_type::shell_id_type const&)>(&impl_type::get_shell))
        .def("num_domains_per_type", &impl_type::num_domains_per_type)
        .def("__len__", &impl_type::num_domains)
        .def("__getitem__", &impl_type::get_domain)
        .def("__iter__", &impl_type::get_domains,
                return_value_policy<return_by_value>())
        ;

    peer::wrappers::generator_wrapper<
        ptr_generator<typename impl_type::domain_id_pair_generator,
                      std::auto_ptr<
                        typename impl_type::domain_id_pair_generator> > >::__register_class("DomainIDPairGenerator");

    peer::converters::register_tuple_converter<typename impl_type::domain_id_pair>();
    peer::converters::register_tuple_converter<get_shell_result_type>();
    shell_variant_converter<Timpl>::__register();
}

} // namespace binding
#endif /* EGFRD_SIMULATOR_HPP */
