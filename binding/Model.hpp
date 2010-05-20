#ifndef BINDING_MODEL_HPP 
#define BINDING_MODEL_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Timpl>
static void Model___setitem__(Timpl* model, std::string const& key, std::string const& value)
{
    (*model)[key] = value;
}

template<typename Tmodel_>
inline void register_model_class(char const* name)
{
    using namespace boost::python;
    typedef Tmodel_ impl_type;

    class_<impl_type, boost::noncopyable>(name)
        .add_property("network_rules",
            make_function(&impl_type::network_rules,
                return_value_policy<reference_existing_object>()))
        .def("add_species_type", &impl_type::add_species_type)
        .def("get_species_type_by_id", &impl_type::get_species_type_by_id)
        .def("__getitem__", (std::string const&(impl_type::*)(std::string const&) const)
                &impl_type::operator[], return_value_policy<copy_const_reference>())
        .def("__setitem__", &Model___setitem__<Tmodel_>)
        .add_property("attributes",
                peer::util::range_from_range<
                    typename impl_type::attributes_range,
                    impl_type, &impl_type::attributes>())
        .add_property("species_types",
                peer::util::range_from_range<
                    typename impl_type::species_type_range,
                    impl_type, &impl_type::get_species_types>());
        ;
}

} // namespace binding

#endif /* BINDING_MODEL_HPP */
