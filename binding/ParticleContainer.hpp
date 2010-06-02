#ifndef BINDING_PARTICLE_CONTAINER_HPP
#define BINDING_PARTICLE_CONTAINER_HPP

#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/override.hpp>
#include <boost/python/wrapper.hpp>
#include "peer/compat.h"
#include "peer/wrappers/generator/pyiterator_generator.hpp"
#include "peer/wrappers/generator/generator_wrapper.hpp"
#include "peer/converters/tuple.hpp"
#include "peer/converters/generator/to_python.hpp"

namespace binding {

template<typename Timpl_>
struct particle_id_pair_generator_converter
{
    typedef Timpl_ native_type;

    struct to_native_converter
    {
        static void* convert(PyObject* ptr)
        {
            boost::python::handle<> iter(
                boost::python::allow_null(PyObject_GetIter(ptr)));
            if (!iter)
            {
                boost::python::throw_error_already_set();
            }
            return new peer::wrappers::pyiterator_generator<
                    typename native_type::result_type>(iter);
        }

        static PyTypeObject const* expected_pytype()
        {
            return &PyBaseObject_Type;
        }
    };

    static void __register()
    {
        peer::wrappers::generator_wrapper<
                ptr_generator<native_type, std::auto_ptr<native_type> > >
                ::__register_class("ParticleIDPairGenerator");
        boost::python::to_python_converter<
                native_type*,
                peer::converters::ptr_generator_to_pyiterator_converter<
                    native_type, std::auto_ptr<native_type> > >();
                    

        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};

template<typename Timpl_>
struct particle_id_pair_and_distance_list_converter
{
    typedef Timpl_ native_type;

    struct to_python_converter
    {
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
                c.set(i, v);
            }

            static reference get(T_ const& c, size_type i)
            {
                return c.at(i);
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

        typedef peer::wrappers::stl_container_wrapper<native_type, std::auto_ptr<native_type>, peer::wrappers::default_policy_generator<policy> > wrapper_type;
        static PyObject* convert(native_type* v)
        {
            return reinterpret_cast<PyObject*>(wrapper_type::create(std::auto_ptr<native_type>(v ? v: new native_type())));
        }
    };

    struct to_native_converter
    {
        static PyTypeObject const* expected_pytype()
        {
            return &PyBaseObject_Type;
        }

        static void* convert(PyObject* pyo)
        {
            if (pyo == Py_None)
            {
                return 0;
            }

            Py_ssize_t hint(-1);
            hint = PyObject_Size(pyo);
            if (hint < 0)
            {
                PyErr_Clear();
            }

            boost::python::handle<> iter(PyObject_GetIter(pyo));
            if (!iter)
            {
                boost::python::throw_error_already_set();
            }

            if (hint < 0)
            {
                hint = compat_PyObject_LengthHint(iter.get());
                if (hint < 0)
                {
                    hint = 0;
                }
            }

            native_type* obj(new native_type);
            obj->reserve(hint);
            while (PyObject* item = PyIter_Next(iter.get()))
            {
                obj->push_back(boost::python::extract<typename native_type::value_type>(item)());
            }

            return obj;
        }
    };

    static void __register()
    {
        to_python_converter::wrapper_type::__class_init__("ParticleIDAndDistanceVector", boost::python::scope().ptr());
        boost::python::to_python_converter<native_type*, to_python_converter>();
        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};

template<typename Tbase_>
class ParticleContainerWrapper
    : public Tbase_, public boost::python::wrapper<Tbase_>
{
public:
    typedef boost::python::wrapper<Tbase_> py_wrapper_type;
    typedef Tbase_ wrapped_type;
    typedef typename wrapped_type::size_type size_type;
    typedef typename wrapped_type::length_type length_type;
    typedef typename wrapped_type::particle_id_type particle_id_type;
    typedef typename wrapped_type::particle_shape_type particle_shape_type;
    typedef typename wrapped_type::position_type position_type;
    typedef typename wrapped_type::species_id_type species_id_type;
    typedef typename wrapped_type::species_type species_type;
    typedef typename wrapped_type::structure_id_type structure_id_type;
    typedef typename wrapped_type::structure_type structure_type;
    typedef typename wrapped_type::particle_id_pair particle_id_pair;
    typedef typename wrapped_type::transaction_type transaction_type;
    typedef typename wrapped_type::particle_id_pair_generator particle_id_pair_generator;
    typedef typename wrapped_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

public:
    virtual ~ParticleContainerWrapper() {}

    virtual size_type num_particles() const
    {
        boost::python::handle<> retval(
            boost::python::allow_null(
                PyObject_GetAttrString(
                    boost::python::detail::wrapper_base_::get_owner(*this),
                    "num_particles")));
        return boost::python::extract<size_type>(retval.get())();
    }

    virtual length_type world_size() const
    {
        boost::python::handle<> retval(
            boost::python::allow_null(
                PyObject_GetAttrString(
                    boost::python::detail::wrapper_base_::get_owner(*this),
                    "world_size")));
        return boost::python::extract<length_type>(retval.get())();
    }

    virtual species_type const& get_species(species_id_type const& id) const
    {
        return py_wrapper_type::get_override("get_species")(id).template unchecked<species_type const&>();
    }

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return py_wrapper_type::get_override("get_structure")(id);
    }

    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        return py_wrapper_type::get_override("new_particle")(sid, pos);
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        return py_wrapper_type::get_override("update_particle")(pi_pair);
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        return py_wrapper_type::get_override("remove_particle")(id);
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        return py_wrapper_type::get_override("get_particle")(id);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_id_pair const& s) const
    {
        return py_wrapper_type::get_override("check_overlap")(
                s.second.shape(), boost::python::make_tuple(s.first))
                .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return py_wrapper_type::get_override("check_overlap")(s, boost::python::tuple())
                .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return py_wrapper_type::get_override("check_overlap")(
                s, boost::python::make_tuple(ignore))
               .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return py_wrapper_type::get_override("check_overlap")(
                s, boost::python::make_tuple(ignore1, ignore2))
               .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_generator* get_particles() const
    {
        return py_wrapper_type::get_override("__iter__")()
                .template unchecked<particle_id_pair_generator*>();
    }

    virtual transaction_type* create_transaction()
    {
        return py_wrapper_type::get_override("create_transaction")()
                .template unchecked<transaction_type*>();
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return py_wrapper_type::get_override("distance")(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return py_wrapper_type::get_override("apply_boundary")(v);
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return py_wrapper_type::get_override("apply_boundary")(v);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return py_wrapper_type::get_override("cyclic_transpose")(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return py_wrapper_type::get_override("cyclic_transpose")(p0, p1);
    }
};


template<typename Timpl>
inline boost::python::objects::class_base register_particle_container_class(
        char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;
    peer::converters::register_tuple_converter<typename impl_type::particle_id_pair>();
    peer::converters::register_tuple_converter<typename impl_type::particle_id_pair_and_distance>();

    particle_id_pair_and_distance_list_converter<typename impl_type::particle_id_pair_and_distance_list>::__register();
    particle_id_pair_generator_converter<typename impl_type::particle_id_pair_generator>::__register();

    return class_<ParticleContainerWrapper<impl_type>, boost::noncopyable>(name)
        .add_property("num_particles", &impl_type::num_particles)
        .add_property("world_size", &impl_type::world_size)
        .def("get_species",
            pure_virtual((typename impl_type::species_type const&(impl_type::*)(typename impl_type::species_id_type const&) const)&impl_type::get_species),
            return_internal_reference<>())
        .def("get_structure", pure_virtual(&impl_type::get_structure))
        .def("new_particle", pure_virtual(&impl_type::new_particle))
        .def("update_particle", pure_virtual(&impl_type::update_particle))
        .def("remove_particle", pure_virtual(&impl_type::remove_particle))
        .def("get_particle", pure_virtual(&impl_type::get_particle))
        .def("check_overlap", pure_virtual((typename impl_type::particle_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&, typename impl_type::particle_id_type const&) const)&impl_type::check_overlap), return_value_policy<return_by_value>())
        .def("check_overlap", pure_virtual((typename impl_type::particle_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&, typename impl_type::particle_id_type const&, typename impl_type::particle_id_type const&) const)&impl_type::check_overlap), return_value_policy<return_by_value>())
        .def("check_overlap", pure_virtual((typename impl_type::particle_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&) const)&impl_type::check_overlap), return_value_policy<return_by_value>())
        .def("create_transaction", pure_virtual(&impl_type::create_transaction),
                return_value_policy<manage_new_object>())
        .def("distance", pure_virtual(&impl_type::distance))
        .def("apply_boundary", pure_virtual((typename impl_type::position_type(impl_type::*)(typename impl_type::position_type const&) const)&impl_type::apply_boundary))
        .def("apply_boundary", pure_virtual((typename impl_type::length_type(impl_type::*)(typename impl_type::length_type const&) const)&impl_type::apply_boundary))
        .def("cyclic_transpose", pure_virtual((typename impl_type::position_type(impl_type::*)(typename impl_type::position_type const&, typename impl_type::position_type const&) const)&impl_type::cyclic_transpose))
        .def("cyclic_transpose", pure_virtual((typename impl_type::length_type(impl_type::*)(typename impl_type::length_type const&, typename impl_type::length_type const&) const)&impl_type::cyclic_transpose))
        .def("__iter__", pure_virtual(&impl_type::get_particles),
                return_value_policy<return_by_value>())
        ;
}

} // namespace binding

#endif /* BINDING_PARTICLE_CONTAINER_HPP */
