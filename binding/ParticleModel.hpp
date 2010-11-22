#ifndef BINDING_MODEL_HPP 
#define BINDING_MODEL_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Tparticle_model_>
inline void register_particle_model_class(char const* name)
{
    using namespace boost::python;
    typedef Tparticle_model_ impl_type;

    class_<impl_type, bases<typename impl_type::base_type>, boost::noncopyable>(name)
        .def("add_structure_type", &impl_type::add_structure_type)
        .def("get_structure_type_by_id", &impl_type::get_structure_type_by_id)
        .add_property("structure_types",
                peer::util::range_from_range<
                    typename impl_type::structure_type_range,
                    impl_type, &impl_type::get_structure_types>())
        ;
}

} // namespace binding

#endif /* BINDING_MODEL_HPP */
