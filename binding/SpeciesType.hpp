#ifndef BINDING_SPECIES_TYPE_HPP
#define BINDING_SPECIES_TYPE_HPP

#include <boost/python.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/lexical_cast.hpp>

namespace binding {

template<typename Timpl_>
struct SpeciesTypeExtras
{
    typedef Timpl_ impl_type;

    static std::string const& __getitem__(impl_type* impl, std::string const& key)
    {
        return (*impl)[key];
    }

    static void __setitem__(impl_type* impl, std::string const& key, std::string const& val)
    {
        (*impl)[key] = val;
    }

    static std::string __str__(impl_type* impl)
    {
        return boost::lexical_cast<std::string>(*impl);
    }
};

template<typename Timpl_>
static boost::python::objects::class_base register_species_type_class(
        char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;
    typedef SpeciesTypeExtras<impl_type> extras_type;

    return class_<impl_type, impl_type*>(name, no_init)
        .add_property("id",
                make_function(&impl_type::id,
                    return_value_policy<copy_const_reference>()))
        .def("__str__", &extras_type::__str__)
        .def("__getitem__", &extras_type::__getitem__,
                return_value_policy<copy_const_reference>())
        .def("__setitem__", &extras_type::__setitem__)
        ;
}

} // namespace binding

#endif /* BINDING_SPECIES_TYPE_HPP */
