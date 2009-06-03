#ifndef PEER_SPECIES_TYPE_HPP
#define PEER_SPECIES_TYPE_HPP

#include <boost/python.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/lexical_cast.hpp>
#include "../SpeciesType.hpp"

namespace peer {

class SpeciesType
{
public:
    typedef ::SpeciesType impl_type;

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

public:
    static void __register_class()
    {
        using namespace boost::python;
        class_<impl_type, impl_type*>("SpeciesType", no_init)
            .add_property("id",
                    make_function(&impl_type::id,
                        return_value_policy<copy_const_reference>()))
            .def("__str__", &SpeciesType::__str__)
            .def("__getitem__", &SpeciesType::__getitem__,
                    return_value_policy<copy_const_reference>())
            .def("__setitem__", &SpeciesType::__setitem__)
            ;
    }
};

} // namespace peer

#endif /* PEER_SPECIES_TYPE_HPP */
