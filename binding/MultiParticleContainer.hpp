#ifndef BINDING_MULTI_PARTICLE_CONTAINER_HPP
#define BINDING_MULTI_PARTICLE_CONTAINER_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl, typename Tworld>
boost::python::objects::class_base register_multi_particle_container_class(char const* name)
{
    using namespace boost::python;

    return class_<Timpl, bases<typename Tworld::particle_container_type>,
           boost::noncopyable>(name, init<Tworld&>());
}

} // namespace binding

#endif /* BINDING_MULTI_PARTICLE_CONTAINER_HPP */
