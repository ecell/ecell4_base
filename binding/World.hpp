#ifndef BINDING_WORLD_HPP
#define BINDING_WORLD_HPP

#include <set>
#include <boost/python.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/const_iterator.hpp>

#include "peer/utils.hpp"
#include "peer/set_indexing_suite.hpp"
#include "utils/range.hpp"
#include "utils/pair.hpp"

namespace binding {

template<typename Timpl_>
struct species_range_converter: public boost::python::default_call_policies
{
    typedef Timpl_ native_type;

    template<typename T_>
    struct policy
    {
        typedef typename boost::range_size<T_>::type size_type;
        typedef typename boost::range_value<T_>::type value_type;
        typedef value_type const& reference;
        typedef value_type const& const_reference;
        typedef typename boost::range_const_iterator<T_>::type iterator;
        typedef typename boost::range_const_iterator<T_>::type const_iterator;

        static size_type size(T_ const& c)
        {
            return ::size(c);
        }

        static void set(T_& c, size_type i, const_reference v)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
        }

        static reference get(T_ const& c, size_type i)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
            throw 0;
        }

        static iterator begin(T_ const& c)
        {
            return boost::begin(c);
        }

        static iterator end(T_ const& c)
        {
            return boost::end(c);
        }
    };

    struct instance_holder
    {
        instance_holder(native_type const& instance,
                        boost::python::handle<> owner)
            : instance_(instance), owner_(owner) {}

        native_type const& operator*() const
        {
            return instance_;
        }

        native_type const* operator->() const
        {
            return &(**this);
        }

        native_type& operator*()
        {
            PyErr_SetString(PyExc_RuntimeError, "object is immutable");
            boost::python::throw_error_already_set();
            return *static_cast<native_type*>(0);
        }

        native_type* operator->()
        {
            return &(**this);
        }

        native_type instance_;
        boost::python::handle<> owner_;
    };

    typedef peer::wrappers::stl_container_wrapper<native_type, instance_holder, peer::wrappers::default_policy_generator<policy> > wrapper_type;

    struct result_converter
    {
        template<typename T_>
        struct apply
        {
            struct type {
                PyObject* operator()(native_type const& val) const
                {
                    return wrapper_type::create(instance_holder(val, boost::python::handle<>()));
                }

                PyTypeObject const* get_pytype() const
                {
                    return &wrapper_type::__class__;
                }
            };
        };
    };

    template<typename Targs>
    static PyObject* postcall(Targs const& arg, PyObject* result)
    {
        reinterpret_cast<wrapper_type*>(result)->ptr().owner_ = boost::python::handle<>(boost::python::borrowed(PyTuple_GET_ITEM(arg, 0)));
        return result;
    }

    static void __register()
    {
        wrapper_type::__class_init__("SpeciesRange", boost::python::scope().ptr());
    }
};

template<typename Timpl_>
struct surfaces_range_converter: public boost::python::default_call_policies
{
    typedef Timpl_ native_type;

    template<typename T_>
    struct policy
    {
        typedef typename boost::range_size<T_>::type size_type;
        typedef typename boost::range_value<T_>::type value_type;
        typedef value_type const& reference;
        typedef value_type const& const_reference;
        typedef typename boost::range_const_iterator<T_>::type iterator;
        typedef typename boost::range_const_iterator<T_>::type const_iterator;

        static size_type size(T_ const& c)
        {
            return ::size(c);
        }

        static void set(T_& c, size_type i, const_reference v)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
        }

        static reference get(T_ const& c, size_type i)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
            throw 0;
        }

        static iterator begin(T_ const& c)
        {
            return boost::begin(c);
        }

        static iterator end(T_ const& c)
        {
            return boost::end(c);
        }
    };

    struct instance_holder
    {
        instance_holder(native_type const& instance,
                        boost::python::handle<> owner)
            : instance_(instance), owner_(owner) {}

        native_type const& operator*() const
        {
            return instance_;
        }

        native_type const* operator->() const
        {
            return &(**this);
        }

        native_type& operator*()
        {
            PyErr_SetString(PyExc_RuntimeError, "object is immutable");
            boost::python::throw_error_already_set();
            return *static_cast<native_type*>(0);
        }

        native_type* operator->()
        {
            return &(**this);
        }

        native_type instance_;
        boost::python::handle<> owner_;
    };

    typedef peer::wrappers::stl_container_wrapper<native_type, instance_holder, peer::wrappers::default_policy_generator<policy> > wrapper_type;

    struct result_converter
    {
        template<typename T_>
        struct apply
        {
            struct type {
                PyObject* operator()(native_type const& val) const
                {
                    return wrapper_type::create(instance_holder(val, boost::python::handle<>()));
                }

                PyTypeObject const* get_pytype() const
                {
                    return &wrapper_type::__class__;
                }
            };
        };
    };

    template<typename Targs>
    static PyObject* postcall(Targs const& arg, PyObject* result)
    {
        reinterpret_cast<wrapper_type*>(result)->ptr().owner_ = boost::python::handle<>(boost::python::borrowed(PyTuple_GET_ITEM(arg, 0)));
        return result;
    }

    static void __register()
    {
        wrapper_type::__class_init__("SurfacesRange", boost::python::scope().ptr());
    }
};

template<typename T>
static typename get_select_first_range<typename T::particle_id_pair_range>::type
World_get_particle_ids(T const& world)
{
    return make_select_first_range(world.get_particles_range());
}


template<typename Timpl_, typename Tbase_, typename Ttraits_>
inline boost::python::objects::class_base register_world_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;
    typedef species_range_converter<typename impl_type::species_range> species_range_converter_type;
    typedef surfaces_range_converter<typename impl_type::surfaces_range> surfaces_range_converter_type;

    species_range_converter_type::__register();
    surfaces_range_converter_type::__register();

    class_<std::set<typename impl_type::particle_id_type> >("ParticleIDSet")
        .def(peer::util::set_indexing_suite<std::set<typename impl_type::particle_id_type> >())
        ;

    return class_<impl_type, bases<Tbase_> >(
        "World", init<typename impl_type::length_type, typename impl_type::size_type>())
        .add_property("cell_size", &impl_type::cell_size)
        .add_property("matrix_size", &impl_type::matrix_size)
        .add_property("species",
            make_function(
                (typename impl_type::species_range(impl_type::*)() const)&impl_type::get_species, species_range_converter_type()))
        .add_property("surfaces",
            make_function(
                (typename impl_type::surfaces_range(impl_type::*)() const)&impl_type::get_surfaces, surfaces_range_converter_type()))
        .add_property("particle_ids", &World_get_particle_ids<impl_type>)
        .def("add_species", &impl_type::add_species)
        .def("add_surface", &impl_type::add_surface)
		.def("get_particle_ids", &impl_type::get_particle_ids)
        .def("distance", &impl_type::template distance<typename impl_type::position_type>)
        .def("distance", &impl_type::template distance<typename Ttraits_::sphere_type>)
        .def("distance", &impl_type::template distance<typename Ttraits_::cylinder_type>)
        .def("distance", &impl_type::template distance<typename Ttraits_::box_type>)
        .def("calculate_pair_CoM", &impl_type::template calculate_pair_CoM<typename impl_type::position_type>)
        ;
}

} // namespace binding

#endif /* BINDING_WORLD_HPP */
