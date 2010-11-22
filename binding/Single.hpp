#ifndef BINDING_SINGLE_HPP
#define BINDING_SINGLE_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
static void single_set_particle(Timpl& impl, typename Timpl::particle_id_pair const& pair)
{
    BOOST_ASSERT(impl.particle().first == pair.first);
    impl.particle().second = pair.second;
}

template<typename Timpl>
inline boost::python::objects::class_base register_single_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .add_property("particle",
            make_function(
                static_cast<typename impl_type::particle_id_pair const&(impl_type::*)() const>(&impl_type::particle),
                return_value_policy<return_by_value>()),
            make_function(&single_set_particle<impl_type>))
        ;
}

} // namespace binding

#endif /* BINDING_SINGLE_HPP */
