#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <exception>
#include <stdexcept>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <boost/lexical_cast.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/refcount.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/converter/object_manager.hpp>

#include <numpy/arrayobject.h>

#include "peer/utils.hpp"
#include "peer/py_hash_support.hpp"
#include "peer/converters/tuple.hpp"
#include "peer/converters/sequence.hpp"
#include "peer/set_indexing_suite.hpp"
#include "peer/numpy/ndarray_converters.hpp"
#include "peer/numpy/wrapped_multi_array.hpp"
#include "peer/numpy/scalar_converters.hpp"
#include "peer/Particle.hpp"
#include "peer/ShellWrapper.hpp"
#include "peer/MatrixSpace.hpp"
#include "peer/SpeciesType.hpp"
#include "peer/Identifier.hpp"
#include "peer/ReactionRule.hpp"
#include "peer/Exception.hpp"
#include "peer/RandomNumberGenerator.hpp"
#include "peer/wrappers/range/pyiterable_range.hpp"
#include "peer/wrappers/range/stl_container_wrapper.hpp"
#include "peer/wrappers/generator/generator_wrapper.hpp"
#include "peer/converters/generator/to_python.hpp"
#include "PyEventScheduler.hpp"

#include "utils/array_traits.hpp"
#include "utils.hpp"
#include "gsl_rng_base.hpp"
#include "geometry.hpp"
#include "MatrixSpace.hpp"
#include "Vector3.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Point.hpp"
#include "Model.hpp"
#include "World.hpp"
#include "GSLRandomNumberGenerator.hpp"
#include "EGFRDSimulator.hpp"
#include "BDPropagator.hpp"
#include "StructureUtils.hpp"
#include "AnalyticalSingle.hpp"
#include "AnalyticalPair.hpp"

typedef CyclicWorldTraits<Real, Real> world_traits_type;
typedef World<world_traits_type> CyclicWorld;
typedef Model model_type;
typedef CyclicWorld::transaction_type transaction_type;
typedef CyclicWorld::base_type::base_type particle_container_type;
typedef EGFRDSimulatorTraitsBase<CyclicWorld> egfrd_simulator_traits_type;
typedef BDPropagator<egfrd_simulator_traits_type> _BDPropagator;

static boost::python::object species_type_class;
static boost::python::object species_info_class;

typedef peer::util::detail::array_to_ndarray_converter<world_traits_type::position_type, world_traits_type::position_type::value_type, 3> position_to_ndarray_converter;

struct ndarray_to_position_converter
{
    typedef world_traits_type::position_type native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PyArray_Check(ptr))
        {
            return NULL;
        }
        
        PyObject* retval(PyArray_CastToType(
                             reinterpret_cast<PyArrayObject*>(ptr),
                             PyArray_DescrFromType(
                                 peer::util::get_numpy_typecode<
                                     native_type::value_type >::value ), 0) );
        if (!retval)
        {
            return NULL;
        }
        
        if (PyObject_Size(retval) != 3)
        {
            boost::python::decref(retval);
            return NULL;
        }
        
        return retval;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        PyArrayObject* array_obj = static_cast<PyArrayObject*>(
            data->stage1.convertible);
        data->stage1.convertible = new(data->storage.bytes) native_type(
            reinterpret_cast<native_type::value_type*>(PyArray_DATA(array_obj)));
        boost::python::decref(reinterpret_cast<PyObject*>(array_obj));
    }
};

struct tuple_to_position_converter
{
    typedef world_traits_type::position_type native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PyTuple_Check(ptr))
        {
            return NULL;
        }
        
        if (PyTuple_GET_SIZE(ptr) != 3)
        {
            return NULL;
        }
        
        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        data->stage1.convertible = new(data->storage.bytes) native_type(
            PyFloat_AsDouble(PyTuple_GET_ITEM(ptr, 0)),
            PyFloat_AsDouble(PyTuple_GET_ITEM(ptr, 1)),
            PyFloat_AsDouble(PyTuple_GET_ITEM(ptr, 2)));
    }
};

struct list_to_position_converter
{
    typedef world_traits_type::position_type native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PyList_Check(ptr))
        {
            return NULL;
        }
        
        if (PyList_GET_SIZE(ptr) != 3)
        {
            return NULL;
        }
        
        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        data->stage1.convertible = new(data->storage.bytes) native_type(
            PyFloat_AsDouble(PyList_GET_ITEM(ptr, 0)),
            PyFloat_AsDouble(PyList_GET_ITEM(ptr, 1)),
            PyFloat_AsDouble(PyList_GET_ITEM(ptr, 2)));
    }
};

struct sphere_to_python_converter
{
    typedef Sphere<world_traits_type::length_type> native_type;

    static PyObject* convert(native_type const& v)
    {
        return boost::python::incref(
            boost::python::object(boost::make_tuple(
                v.position(), v.radius())).ptr());
    }
};

struct python_to_sphere_converter
{
    typedef Sphere<world_traits_type::length_type> native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyTuple_Check(pyo))
        {
            return 0;
        }
        if (PyTuple_Size(pyo) != 2)
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        PyObject* items[] = { PyTuple_GetItem(pyo, 0), PyTuple_GetItem(pyo, 1) };
        void* storage(reinterpret_cast<
            boost::python::converter::rvalue_from_python_storage<
                native_type >*
            >(data)->storage.bytes);
        new (storage) native_type (
            boost::python::extract<
                native_type::position_type>(items[0]),
            PyFloat_AsDouble(items[1]));
        data->convertible = storage;
    }
};

struct species_type_to_species_id_converter
{
    typedef world_traits_type::species_id_type native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, reinterpret_cast<PyTypeObject*>(
                species_type_class.ptr())))
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;
        void* storage(reinterpret_cast<
            converter::rvalue_from_python_storage<native_type>* >(
                data)->storage.bytes);
        new (storage) native_type(static_cast<SpeciesType*>(extract<SpeciesType*>(object(borrowed(pyo))))->id());
        data->convertible = storage;
    }
};

struct species_info_to_species_id_converter
{
    typedef world_traits_type::species_type::identifier_type native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, reinterpret_cast<PyTypeObject*>(
                species_info_class.ptr())))
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;
        void* storage(reinterpret_cast<
            converter::rvalue_from_python_storage<native_type>* >(
                data)->storage.bytes);
        new (storage) native_type(static_cast<world_traits_type::species_type*>(extract<world_traits_type::species_type*>(object(borrowed(pyo))))->id());
        data->convertible = storage;
    }
};

static void Model___setitem__(Model* model, std::string const& key, std::string const& value)
{
    (*model)[key] = value;
}


template<typename T_>
static void register_id_generator(char const* class_name)
{
    using namespace boost::python;

    class_<SerialIDGenerator<T_> >(class_name, init<int>())
        .def("__call__", &SerialIDGenerator<T_>::operator())
        ;
}

template<typename T_>
static T_ 
calculate_pair_CoM(T_ const& p1, 
                   T_ const& p2, 
                   typename element_type_of< T_ >::type const& D1,
                   typename element_type_of< T_ >::type const& D2,
                   typename element_type_of< T_ >::type const& world_size)
{
    typedef typename element_type_of<T_>::type element_type;   

    T_ retval;

    const T_ p2t(cyclic_transpose<T_>(p2, p1, world_size));

    return modulo(
        divide(
            add(multiply(p1, D2), multiply(p2t, D1)),
            add(D1, D2)),
        world_size);
}

template<typename T_>
static boost::python::object Sphere___getitem__(T_ const& obj, int index)
{
    switch (index)
    {
    default:
        PyErr_SetString(PyExc_IndexError, "index out of range");
        boost::python::throw_error_already_set();
        break;
    case 0:
        return boost::python::object(obj.position());
    case 1:
        return boost::python::object(obj.radius());
    }

    return boost::python::object();
}

struct reaction_rule_vector_converter
{
    typedef egfrd_simulator_traits_type::network_rules_type::reaction_rule_vector native_type;

    struct instance_holder
    {
        instance_holder(native_type const& instance): instance_(instance) {}

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

        native_type const& instance_;
    };

    typedef peer::wrappers::stl_container_wrapper<native_type, instance_holder> wrapper_type;

    struct to_python_converter
    {
        static PyObject* convert(native_type const& v)
        {
            return reinterpret_cast<PyObject*>(wrapper_type::create(instance_holder(v)));
        }
    };

    struct to_native_converter
    {
        static void* convert(PyObject* ptr)
        {
            return const_cast<native_type*>(&*reinterpret_cast<wrapper_type const*>(ptr)->ptr());
        }

        static PyTypeObject const* expected_pytype()
        {
            return &wrapper_type::__class__;
        }
    };

    static void __register()
    {
        wrapper_type::__class_init__("ReactionRuleVector", boost::python::scope().ptr());
        boost::python::to_python_converter<native_type, to_python_converter>();
        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};

struct particle_id_pair_and_distance_list_converter
{
    typedef CyclicWorld::particle_id_pair_and_distance_list native_type;

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
                hint = _PyObject_LengthHint(iter.get(), 0);
                if (hint < 0)
                {
                    hint = 0;
                }
            }

            native_type* obj(new native_type);
            obj->reserve(hint);
            while (PyObject* item = PyIter_Next(iter.get()))
            {
                obj->push_back(boost::python::extract<native_type::value_type>(item)());
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

struct particle_id_pair_range_select_first_converter
{
    typedef get_select_first_range<CyclicWorld::particle_id_pair_range>::type native_type;

    struct holder
    {
        holder(native_type const& val): val_(val) {}

        operator native_type const*() const
        {
            return &val_;
        }

        operator native_type*() const
        {
            // this should never be called actually
            return &const_cast<holder*>(this)->val_;
        }

        native_type const& operator*() const
        {
            return val_;
        }

    private:
        native_type val_;
    };

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

        static void set(T_ const& c, size_type i, const_reference v)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
        }

        static value_type get(T_ const& c, size_type i)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
            throw 0; // dummy
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

    typedef peer::wrappers::stl_container_wrapper<native_type, holder, peer::wrappers::default_policy_generator<policy> > wrapper_type;

    struct to_python_converter
    {
        static PyObject* convert(native_type const& v)
        {
            return reinterpret_cast<PyObject*>(wrapper_type::create(v));
        }
    };

    struct to_native_converter
    {
        static void* convertible(PyObject* ptr)
        {
            if (Py_TYPE(ptr) != &wrapper_type::__class__)
            {
                return NULL;
            }
            return ptr;
        }
        
        static void construct(PyObject* ptr,
                              boost::python::converter::rvalue_from_python_storage<native_type>* data)
        {
            data->stage1.convertible = new(data->storage.bytes) native_type(
                *reinterpret_cast<wrapper_type*>(data->stage1.convertible)->ptr());
        }
    };

    static void __register()
    {
        wrapper_type::__class_init__("ParticleIDRange", boost::python::scope().ptr());
        boost::python::to_python_converter<native_type, to_python_converter>();
        peer::util::to_native_converter<native_type, to_native_converter>();
    }
};

struct species_range_converter: public boost::python::default_call_policies
{
    typedef CyclicWorld::species_range native_type;

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

struct surfaces_range_converter: public boost::python::default_call_policies
{
    typedef CyclicWorld::surfaces_range native_type;

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

struct domain_id_pair_converter
{
    typedef egfrd_simulator_traits_type::domain_id_pair native_type;

    struct to_python_converter
    {
        static PyObject* convert(native_type const& v)
        {
            return boost::python::incref(boost::python::make_tuple(
                v.first, v.second).ptr());
        }
    };

    static void __register()
    {
        boost::python::to_python_converter<native_type, to_python_converter>();
        peer::util::to_native_converter<native_type,
            peer::converters::pytuple_to_tuple_converter<native_type> >();
    }
};

struct particle_id_pair_generator_converter
{
    typedef CyclicWorld::particle_id_pair_generator native_type;

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
                    native_type::result_type>(iter);
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

template<typename T>
static typename T::particle_id_pair_and_distance_list* World_check_overlap(
    T const& world,
    typename T::particle_shape_type const& s,
    twofold_container<typename T::particle_id_type> const& ignore)
{
    return world.check_overlap(s, ignore);
}

template<typename T>
static typename get_select_first_range<typename T::particle_id_pair_range>::type
World_get_particle_ids(T const& world)
{
    return make_select_first_range(world.get_particles_range());
}

struct static_gsl_rng: public gsl_rng_base<static_gsl_rng>
{
    typedef unsigned int result_type;

    static const char name[];
    static const unsigned long int min = 0;
    static const unsigned long int max = static_cast<result_type>(-1);

    void set(unsigned long int seed)
    {
        idx_ = std::min(seed, static_cast<unsigned long int>(PyObject_Size(seq_.ptr())));
    }

    unsigned long int get()
    {
        Py_ssize_t nelems(PyObject_Size(seq_.ptr()));
        if (idx_ >= nelems)
        {
            return min; 
        }
        boost::python::handle<> i(
                boost::python::allow_null(
                    PySequence_GetItem(seq_.ptr(), idx_)));
        if (!i)
        {
            return min;
        }
        ++idx_;
        return PyLong_AsUnsignedLong(i.get());
    }

    double get_double()
    {
        return static_cast<double>(get()) / (static_cast<double>(max) + 1);
    }

    static_gsl_rng(boost::python::object seq)
        : seq_(seq), idx_(0) {}

private:
    boost::python::object seq_;
    Py_ssize_t idx_;
};

const char static_gsl_rng::name[] = "static_gsl_rng";

template<gsl_rng_type const*& Prng_>
static GSLRandomNumberGenerator create_gsl_rng()
{
    return GSLRandomNumberGenerator(gsl_rng_alloc(Prng_));
}

static GSLRandomNumberGenerator create_static_gsl_rng(boost::python::object seq)
{
    return GSLRandomNumberGenerator(
            GSLRandomNumberGenerator::rng_handle(new static_gsl_rng(seq)));
}

world_traits_type::position_type draw_bd_displacement(world_traits_type::surface_type const& surface, world_traits_type::length_type const& length, egfrd_simulator_traits_type::rng_type& rng)
{
    return StructureUtils<egfrd_simulator_traits_type>::draw_bd_displacement(surface, length, rng);
}

class _EGFRDSimulator: public EGFRDSimulator<egfrd_simulator_traits_type>, public boost::python::wrapper<EGFRDSimulator<egfrd_simulator_traits_type> >
{
    typedef EGFRDSimulator<egfrd_simulator_traits_type> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::network_rules_type network_rules_type;
    typedef base_type::rng_type rng_type;
    typedef base_type::domain_id_pair_generator domain_id_pair_generator;

public:
    _EGFRDSimulator(world_type& world, rng_type& rng,
                    network_rules_type const& network_rules)
        : base_type(world, rng, network_rules) {}

    virtual void step()
    {
        get_override("step")();
    }
};

class _ParticleContainer: public particle_container_type,
                          public boost::python::wrapper<particle_container_type>
{
public:
    virtual ~_ParticleContainer() {}

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
        return get_override("get_species")(id).unchecked<species_type const&>();
    }

    virtual boost::shared_ptr<surface_type> get_surface(surface_id_type const& id) const
    {
        return get_override("get_surface")(id);
    }

    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        return get_override("new_particle")(sid, pos);
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        return get_override("update_particle")(pi_pair);
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        return get_override("remove_particle")(id);
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        return get_override("get_particle")(id);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_id_pair const& s) const
    {
        return get_override("check_overlap")(
                s.second.shape(), boost::python::make_tuple(s.first))
                .unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return get_override("check_overlap")(s, boost::python::tuple())
                .unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return get_override("check_overlap")(
                s, boost::python::make_tuple(ignore))
               .unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return get_override("check_overlap")(
                s, boost::python::make_tuple(ignore1, ignore2))
               .unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_generator* get_particles() const
    {
        return get_override("__iter__")().unchecked<particle_id_pair_generator*>();
    }

    virtual transaction_type* create_transaction()
    {
        return get_override("create_transaction")().unchecked<transaction_type*>();
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return get_override("distance")(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return get_override("apply_boundary")(v);
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return get_override("apply_boundary")(v);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return get_override("cyclic_transpose")(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return get_override("cyclic_transpose")(p0, p1);
    }
};

static void _BDPropagator_propagate_all(_BDPropagator& self)
{
    while (self());
}

static particle_container_type::particle_id_pair_and_distance_list* do_check_overlap(particle_container_type const& pc, particle_container_type::particle_type::shape_type const& x)
{
    return pc.check_overlap(x);
}

BOOST_PYTHON_MODULE( _gfrd )
{
    using namespace boost::python;

    import_array();

    // GSL error handler: is this the best place for this?
    gsl_set_error_handler( &gsl_error_handler );

  
    peer::util::register_std_exception_translator();

    peer::util::register_seq_wrapped_multi_array_converter<world_traits_type::length_type>();
    // peer::util::register_ndarray_wrapped_multi_array_converter<length_type, 1>();
    peer::util::register_ndarray_wrapped_multi_array_converter<world_traits_type::length_type, 2>();
    peer::util::register_ndarray_wrapped_multi_array_converter<world_traits_type::length_type, 3>();

    //peer::util::register_tuple_converter< boost::tuple< Real, EventType > >();

    class_<PyEvent>( "PyEvent", 
                     init<const Real, const object&, const object&>() )
        .def( "setTime", &PyEvent::setTime )
        .def( "getTime", &PyEvent::getTime )
        .def( "getObj", &PyEvent::getObj,
              return_value_policy<copy_const_reference>() )
        .def( "getArg", &PyEvent::getArg,
              return_value_policy<copy_const_reference>() )
        ;


    typedef const PyEventScheduler::Event& 
        (PyEventScheduler::*geteventrefsig)() const;
    typedef const PyEventScheduler::Event& 
        (PyEventScheduler::*geteventrefbyindexsig)
        ( const PyEventScheduler::EventID );

    class_<PyEventScheduler, boost::noncopyable>( "EventScheduler" )
        .def( "getTime", &PyEventScheduler::getTime )
        .def( "getTopTime", &PyEventScheduler::getTopTime )
        .def( "getSize", &PyEventScheduler::getSize )
        .def( "getTopEvent", geteventrefsig( &PyEventScheduler::getTopEvent ),
              return_value_policy<copy_const_reference>() )
        .def( "getTopID", &PyEventScheduler::getTopID )
        .def( "peekSecondEvent", 
              geteventrefsig( &PyEventScheduler::peekSecondEvent ),
              return_value_policy<copy_const_reference>() )
        .def( "getEvent", geteventrefbyindexsig( &PyEventScheduler::getEvent ),
              return_value_policy<copy_const_reference>() )
        .def( "getEventByIndex", &PyEventScheduler::getEventByIndex,
              return_value_policy<copy_const_reference>() )
        .def( "step", &PyEventScheduler::step )
        .def( "clear", &PyEventScheduler::clear )
        .def( "addEvent", &PyEventScheduler::addEvent )
        .def( "removeEvent", &PyEventScheduler::removeEvent )
        .def( "updateEventTime", &PyEventScheduler::updateEventTime )
//        .def( "updateAllEventDependency", 
//              &PyEventScheduler::updateAllEventDependency )
        .def( "check", &PyEventScheduler::check )
        ;


    def( "length_sq", &length_sq< world_traits_type::position_type > );
    def( "length", &length< world_traits_type::position_type > );
    def( "distance", (world_traits_type::length_type(*)(world_traits_type::position_type const&, world_traits_type::position_type const&))&distance<world_traits_type::position_type> );
    def( "distance_cyclic", &distance_cyclic<world_traits_type::position_type, world_traits_type::position_type> );
    def( "apply_boundary", &apply_boundary<world_traits_type::position_type, world_traits_type::length_type> );
    def( "calculate_pair_CoM", &calculate_pair_CoM<world_traits_type::position_type> );

    def( "normalize", (world_traits_type::position_type(*)(world_traits_type::position_type const&))&normalize<world_traits_type::position_type> );
    def( "normalize", (world_traits_type::position_type(*)(world_traits_type::position_type const&, world_traits_type::length_type const&))&normalize<world_traits_type::position_type> );
    def( "cyclic_transpose", &cyclic_transpose<world_traits_type::position_type, element_type_of<world_traits_type::position_type>::type> );

    to_python_converter<world_traits_type::position_type,
        position_to_ndarray_converter>();
#if 0
    to_python_converter<Sphere<world_traits_type::length_type>,
        sphere_to_python_converter>();
#endif
    peer::util::to_native_converter<Sphere<world_traits_type::length_type>,
        python_to_sphere_converter>();
    peer::util::to_native_converter<world_traits_type::position_type,
        ndarray_to_position_converter>();
    peer::util::to_native_converter<world_traits_type::position_type,
        tuple_to_position_converter>();
    peer::util::to_native_converter<world_traits_type::position_type,
        list_to_position_converter>();

    peer::MatrixSpace< MatrixSpace<egfrd_simulator_traits_type::spherical_shell_type, egfrd_simulator_traits_type::shell_id_type> >::__register_class("SphericalShellContainer");
    peer::MatrixSpace< MatrixSpace<egfrd_simulator_traits_type::cylindrical_shell_type, egfrd_simulator_traits_type::shell_id_type> >::__register_class("CylindricalShellContainer");
    peer::MatrixSpace<MatrixSpace<
        world_traits_type::particle_type,
        world_traits_type::particle_id_type, get_mapper_mf> >::__register_class("ParticleContainer");
    species_type_class = peer::SpeciesType::__register_class();

    class_<std::set<world_traits_type::particle_id_type> >("ParticleIDSet")
        .def(peer::util::set_indexing_suite<std::set<world_traits_type::particle_id_type> >())
        ;

    class_<model_type, boost::noncopyable>("Model")
        .add_property("network_rules",
            make_function(&model_type::network_rules,
                return_value_policy<reference_existing_object>()))
        .def("new_species_type", &model_type::new_species_type,
                return_value_policy<reference_existing_object>())
        .def("get_species_type_by_id", &model_type::get_species_type_by_id,
                return_value_policy<reference_existing_object>())
        .def("__getitem__", (std::string const&(model_type::*)(std::string const&) const)
                &model_type::operator[], return_value_policy<copy_const_reference>())
        .def("__setitem__", &Model___setitem__)
        .add_property("attributes",
                peer::util::range_from_range<
                    model_type::attributes_range,
                    model_type, &model_type::attributes>())
        .add_property("species_types",
                peer::util::range_from_range<
                    model_type::species_type_range,
                    model_type, &model_type::get_species_types,
                    reference_existing_object>());
        ;

    peer::ReactionRule::__register_class();

    peer::IdentifierWrapper<world_traits_type::species_id_type>::__register_class("SpeciesTypeID");
    peer::util::to_native_converter<world_traits_type::species_id_type, species_type_to_species_id_converter>();

    peer::wrappers::generator_wrapper<ptr_generator<NetworkRules::reaction_rule_generator, std::auto_ptr<NetworkRules::reaction_rule_generator> > >::__register_class("ReactionRuleGenerator");
    boost::python::to_python_converter<
            NetworkRules::reaction_rule_generator*,
            peer::converters::ptr_generator_to_pyiterator_converter<
                NetworkRules::reaction_rule_generator,
                std::auto_ptr<NetworkRules::reaction_rule_generator> > >();

    peer::util::ExceptionWrapper<not_found, peer::util::PyExcTraits<&PyExc_LookupError> >::__register_class("NotFound");
    peer::util::ExceptionWrapper<already_exists, peer::util::PyExcTraits<&PyExc_StandardError> >::__register_class("AlreadyExists");

    peer::IdentifierWrapper<world_traits_type::particle_id_type>::__register_class("ParticleID");
    register_id_generator<world_traits_type::particle_id_type>("ParticleIDGenerator");
    peer::ParticleWrapper<world_traits_type::particle_type>::__register_class("Particle");

    peer::IdentifierWrapper<egfrd_simulator_traits_type::shell_id_type>::__register_class("ShellID");
    register_id_generator<egfrd_simulator_traits_type::shell_id_type>("ShellIDGenerator");
    peer::ShellWrapper<egfrd_simulator_traits_type::spherical_shell_type>::__register_class("SphericalShell");
    peer::ShellWrapper<egfrd_simulator_traits_type::cylindrical_shell_type>::__register_class("CylindricalShell");

    peer::IdentifierWrapper<egfrd_simulator_traits_type::domain_id_type>::__register_class("DomainID");
    register_id_generator<egfrd_simulator_traits_type::domain_id_type>("DomainIDGenerator");

    class_<NetworkRules, boost::noncopyable>("NetworkRules", no_init)
        .def("add_reaction_rule", &NetworkRules::add_reaction_rule)
        .def("remove_reaction_rule", &NetworkRules::remove_reaction_rule)
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(world_traits_type::species_id_type const&) const>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(world_traits_type::species_id_type const&, world_traits_type::species_id_type const&) const>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        ;


    typedef egfrd_simulator_traits_type::reaction_rule_type reaction_rule_info_type;
    class_<reaction_rule_info_type>("ReactionRuleInfo",
        init<reaction_rule_info_type::identifier_type,
             reaction_rule_info_type::rate_type,
             twofold_container<world_traits_type::species_id_type>,
             std::vector<world_traits_type::species_id_type> >())
        .add_property("id", 
            make_function(&reaction_rule_info_type::id,
                          return_value_policy<return_by_value>()))
        .add_property("k", make_function(&reaction_rule_info_type::k))
        .add_property("products",
            make_function(&reaction_rule_info_type::get_products,
                          return_value_policy<return_by_value>()))
        .add_property("reactants",
            make_function(&reaction_rule_info_type::get_reactants,
                          return_value_policy<return_by_value>()));

    peer::converters::register_range_to_tuple_converter<reaction_rule_info_type::species_id_range>();

    peer::converters::register_iterable_to_range_converter<reaction_rule_info_type::species_id_range>();

    typedef egfrd_simulator_traits_type::network_rules_type NetworkRulesWrapper;
    class_<NetworkRulesWrapper, boost::noncopyable>("NetworkRulesWrapper", init<NetworkRules const&>())
        .def("query_reaction_rule", (NetworkRulesWrapper::reaction_rule_vector const&(NetworkRulesWrapper::*)(NetworkRulesWrapper::species_id_type const&) const)&NetworkRulesWrapper::query_reaction_rule,
            return_value_policy<return_by_value>())
        .def("query_reaction_rule", (NetworkRulesWrapper::reaction_rule_vector const&(NetworkRulesWrapper::*)(NetworkRulesWrapper::species_id_type const&, NetworkRulesWrapper::species_id_type const&) const)&NetworkRulesWrapper::query_reaction_rule,
            return_value_policy<return_by_value>())
        ;

    reaction_rule_vector_converter::wrapper_type::__class_init__("NetworkRulesWrapper.ReactionRuleVector", scope().ptr());
    reaction_rule_vector_converter::__register();

    peer::converters::register_tuple_converter<CyclicWorld::particle_id_pair>();
    peer::converters::register_tuple_converter<CyclicWorld::particle_id_pair_and_distance>();

    particle_id_pair_generator_converter::__register();

    class_<_ParticleContainer, boost::noncopyable>("_ParticleContainer")
        .add_property("num_particles", &particle_container_type::num_particles)
        .add_property("world_size", &particle_container_type::world_size)
        .def("get_species",
            pure_virtual((particle_container_type::species_type const&(particle_container_type::*)(particle_container_type::species_id_type const&) const)&particle_container_type::get_species),
            return_internal_reference<>())
        .def("get_surface", pure_virtual(&particle_container_type::get_surface))
        .def("new_particle", pure_virtual(&particle_container_type::new_particle))
        .def("update_particle", pure_virtual(&particle_container_type::update_particle))
        .def("remove_particle", pure_virtual(&particle_container_type::remove_particle))
        .def("get_particle", pure_virtual(&particle_container_type::get_particle))
        .def("check_overlap", pure_virtual((particle_container_type::particle_id_pair_and_distance_list*(particle_container_type::*)(particle_container_type::particle_type::shape_type const&, particle_container_type::particle_id_type const&) const)&particle_container_type::check_overlap), return_value_policy<return_by_value>())
        .def("check_overlap", pure_virtual((particle_container_type::particle_id_pair_and_distance_list*(particle_container_type::*)(particle_container_type::particle_type::shape_type const&, particle_container_type::particle_id_type const&, particle_container_type::particle_id_type const&) const)&particle_container_type::check_overlap), return_value_policy<return_by_value>())
        .def("check_overlap", pure_virtual((particle_container_type::particle_id_pair_and_distance_list*(particle_container_type::*)(particle_container_type::particle_type::shape_type const&) const)&particle_container_type::check_overlap), return_value_policy<return_by_value>())
        .def("create_transaction", pure_virtual(&particle_container_type::create_transaction),
                return_value_policy<manage_new_object>())
        .def("distance", pure_virtual(&particle_container_type::distance))
        .def("apply_boundary", pure_virtual((particle_container_type::position_type(particle_container_type::*)(particle_container_type::position_type const&) const)&particle_container_type::apply_boundary))
        .def("apply_boundary", pure_virtual((particle_container_type::length_type(particle_container_type::*)(particle_container_type::length_type const&) const)&particle_container_type::apply_boundary))
        .def("cyclic_transpose", pure_virtual((particle_container_type::position_type(particle_container_type::*)(particle_container_type::position_type const&, particle_container_type::position_type const&) const)&particle_container_type::cyclic_transpose))
        .def("cyclic_transpose", pure_virtual((particle_container_type::length_type(particle_container_type::*)(particle_container_type::length_type const&, particle_container_type::length_type const&) const)&particle_container_type::cyclic_transpose))
        .def("__iter__", pure_virtual(&particle_container_type::get_particles),
                return_value_policy<return_by_value>())
        ;

    class_<transaction_type, bases<particle_container_type>, boost::noncopyable>("Transaction", no_init)
        .add_property("added_particles",
            make_function(&transaction_type::get_added_particles,
                return_value_policy<return_by_value>()))
        .add_property("removed_particles",
            make_function(&transaction_type::get_removed_particles,
                return_value_policy<return_by_value>()))
        .add_property("modified_particles",
            make_function(&transaction_type::get_modified_particles,
                return_value_policy<return_by_value>()))
        .def("rollback", &transaction_type::rollback)
        ;

    class_<TransactionImpl<particle_container_type>, bases<transaction_type>, boost::noncopyable>(
           "TransactionImpl", init<particle_container_type&>());

    particle_id_pair_and_distance_list_converter::__register();
    species_range_converter::__register();
    peer::converters::register_range_to_tuple_converter<twofold_container<CyclicWorld::particle_id_type> >();
    peer::converters::register_iterable_to_range_converter<twofold_container<CyclicWorld::particle_id_type> >();
    particle_id_pair_range_select_first_converter::__register();

    class_<CyclicWorld, bases<particle_container_type> >(
        "World", init<CyclicWorld::length_type, CyclicWorld::size_type>())
        .add_property("cell_size", &CyclicWorld::cell_size)
        .add_property("matrix_size", &CyclicWorld::matrix_size)
        .add_property("species",
            make_function(
                (CyclicWorld::species_range(CyclicWorld::*)() const)&CyclicWorld::get_species, species_range_converter()))
        .add_property("surfaces",
            make_function(
                (CyclicWorld::surfaces_range(CyclicWorld::*)() const)&CyclicWorld::get_surfaces, surfaces_range_converter()))
        .add_property("particle_ids", &World_get_particle_ids<CyclicWorld>)
        .def("add_species", &CyclicWorld::add_species)
        .def("add_surface", &CyclicWorld::add_surface)
		.def("get_particle_ids", &CyclicWorld::get_particle_ids)
        .def("distance", &CyclicWorld::distance<CyclicWorld::position_type>)
        .def("distance", &CyclicWorld::distance<egfrd_simulator_traits_type::sphere_type>)
        .def("distance", &CyclicWorld::distance<egfrd_simulator_traits_type::cylinder_type>)
        .def("distance", &CyclicWorld::distance<egfrd_simulator_traits_type::box_type>)
        .def("calculate_pair_CoM", &CyclicWorld::calculate_pair_CoM<CyclicWorld::position_type>)
        ;

    typedef world_traits_type::species_type species_type;

    species_info_class = class_<species_type>("SpeciesInfo",
            init<species_type::identifier_type>())
        .def(init<species_type::identifier_type, species_type::length_type, species_type::D_type, species_type::surface_id_type>())
        .add_property("id",
            make_function(&species_type::id,
                return_value_policy<return_by_value>()))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    species_type, species_type::length_type,
                    &species_type::radius,
                    &species_type::radius>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                species_type, species_type::length_type,
                &species_type::radius,
                &species_type::radius>::set)
        .add_property("surface_id",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    species_type, species_type::surface_id_type,
                    &species_type::surface_id,
                    &species_type::surface_id>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                species_type, species_type::surface_id_type,
                &species_type::surface_id,
                &species_type::surface_id>::set)
        .add_property("D",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    species_type, species_type::D_type,
                    &species_type::D,
                    &species_type::D>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                species_type, species_type::D_type,
                &species_type::D,
                &species_type::D>::set)
        ;

    class_<egfrd_simulator_traits_type::sphere_type>("Sphere")
        .def(init<egfrd_simulator_traits_type::sphere_type::position_type, 
                  egfrd_simulator_traits_type::sphere_type::length_type>())
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::position_type,
                    &egfrd_simulator_traits_type::sphere_type::position,
                    &egfrd_simulator_traits_type::sphere_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::position_type,
                    &egfrd_simulator_traits_type::sphere_type::position,
                    &egfrd_simulator_traits_type::sphere_type::position>::set))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::length_type,
                    &egfrd_simulator_traits_type::sphere_type::radius,
                    &egfrd_simulator_traits_type::sphere_type::radius>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::length_type,
                    &egfrd_simulator_traits_type::sphere_type::radius,
                    &egfrd_simulator_traits_type::sphere_type::radius>::set))
        .def("__getitem__", &Sphere___getitem__<egfrd_simulator_traits_type::sphere_type>);

    class_<egfrd_simulator_traits_type::cylinder_type>("Cylinder")
        .def(init<egfrd_simulator_traits_type::cylinder_type::position_type, 
                  egfrd_simulator_traits_type::cylinder_type::length_type,
                  egfrd_simulator_traits_type::cylinder_type::position_type, 
                  egfrd_simulator_traits_type::cylinder_type::length_type>())
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::position,
                    &egfrd_simulator_traits_type::cylinder_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::position,
                    &egfrd_simulator_traits_type::cylinder_type::position>::set))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::radius,
                    &egfrd_simulator_traits_type::cylinder_type::radius>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::radius,
                    &egfrd_simulator_traits_type::cylinder_type::radius>::set))
        .add_property("size",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::size,
                    &egfrd_simulator_traits_type::cylinder_type::size>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::size,
                    &egfrd_simulator_traits_type::cylinder_type::size>::set))
        .add_property("unit_z",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z>::set));

    to_python_converter<boost::array<egfrd_simulator_traits_type::box_type::length_type, 3>, peer::util::detail::to_ndarray_converter<boost::array<egfrd_simulator_traits_type::box_type::length_type, 3> > >();

    class_<egfrd_simulator_traits_type::box_type>("Box")
        .def(init<egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::length_type,
                  egfrd_simulator_traits_type::box_type::length_type, 
                  egfrd_simulator_traits_type::box_type::length_type>())
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::position,
                    &egfrd_simulator_traits_type::box_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::position,
                    &egfrd_simulator_traits_type::box_type::position>::set))
        .add_property("extent",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    boost::array<egfrd_simulator_traits_type::box_type::length_type, 3>,
                    &egfrd_simulator_traits_type::box_type::extent,
                    &egfrd_simulator_traits_type::box_type::extent>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    boost::array<egfrd_simulator_traits_type::box_type::length_type, 3>,
                    &egfrd_simulator_traits_type::box_type::extent,
                    &egfrd_simulator_traits_type::box_type::extent>::set))
        .add_property("unit_x",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_x,
                    &egfrd_simulator_traits_type::box_type::unit_x>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_x,
                    &egfrd_simulator_traits_type::box_type::unit_x>::set))
        .add_property("unit_y",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_y,
                    &egfrd_simulator_traits_type::box_type::unit_y>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_y,
                    &egfrd_simulator_traits_type::box_type::unit_y>::set))
        .add_property("unit_z",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_z,
                    &egfrd_simulator_traits_type::box_type::unit_z>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_z,
                    &egfrd_simulator_traits_type::box_type::unit_z>::set));

    class_<world_traits_type::surface_type>("Structure", no_init)
        .add_property("id", 
            make_function(&world_traits_type::surface_type::id,
                          return_value_policy<return_by_value>()))
        ;

    typedef egfrd_simulator_traits_type::planar_surface_type PlanarSurface;
    class_<PlanarSurface,
           bases<world_traits_type::surface_type>,
           boost::shared_ptr<PlanarSurface> >(
           "_PlanarSurface", init<PlanarSurface::identifier_type, PlanarSurface::shape_type>())
        .add_property("shape",
            make_function((PlanarSurface::shape_type const&(PlanarSurface::*)()const)&PlanarSurface::shape,
                          return_value_policy<return_by_value>()))
        ;

    typedef egfrd_simulator_traits_type::spherical_surface_type SphericalSurface;
    class_<SphericalSurface,
           bases<world_traits_type::surface_type>,
           boost::shared_ptr<SphericalSurface> >("_SphericalSurface", init<SphericalSurface::identifier_type, SphericalSurface::shape_type>())
        .add_property("shape",
            make_function((SphericalSurface::shape_type const&(SphericalSurface::*)()const)&SphericalSurface::shape,
                          return_value_policy<return_by_value>()))
        ;

    typedef egfrd_simulator_traits_type::cylindrical_surface_type CylindricalSurface;
    class_<CylindricalSurface,
           bases<world_traits_type::surface_type>,
           boost::shared_ptr<CylindricalSurface> >("_CylindricalSurface", init<CylindricalSurface::identifier_type, CylindricalSurface::shape_type>())
        .add_property("shape",
            make_function((CylindricalSurface::shape_type const&(CylindricalSurface::*)()const)&CylindricalSurface::shape,
                          return_value_policy<return_by_value>()))
        ;

    typedef egfrd_simulator_traits_type::cuboidal_region_type CuboidalRegion;
    class_<CuboidalRegion,
           bases<world_traits_type::surface_type>,
           boost::shared_ptr<CuboidalRegion> >("_CuboidalRegion", init<CuboidalRegion::identifier_type, CuboidalRegion::shape_type>())
        .add_property("shape",
            make_function((CuboidalRegion::shape_type const&(CuboidalRegion::*)()const)&CuboidalRegion::shape,
                          return_value_policy<return_by_value>()))
        ;

    typedef Domain<egfrd_simulator_traits_type> Domain;
    class_<Domain>("_Domain", no_init)
        .add_property("surface_id", 
            make_function(&Domain::surface_id,
                return_value_policy<return_by_value>()))
        .add_property("event_id", 
            make_function(
                &peer::util::reference_accessor_wrapper<
                    Domain, Domain::event_id_type,
                    &Domain::event_id, &Domain::event_id>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    Domain, Domain::event_id_type,
                    &Domain::event_id, &Domain::event_id>::set))
        .add_property("last_time", 
            make_function(
                &peer::util::reference_accessor_wrapper<
                    Domain, Domain::time_type,
                    &Domain::last_time, &Domain::last_time>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    Domain, Domain::time_type,
                    &Domain::last_time, &Domain::last_time>::set))
        .add_property("event_kind", 
            make_function(
                &peer::util::reference_accessor_wrapper<
                    Domain, Domain::event_kind_type,
                    &Domain::event_kind, &Domain::event_kind>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    Domain, Domain::event_kind_type,
                    &Domain::event_kind, &Domain::event_kind>::set))
        ;

    typedef _EGFRDSimulator::spherical_single_type SphericalSingle;
    class_<SphericalSingle, bases<Domain>,
           boost::shared_ptr<SphericalSingle> >(
        "_SphericalSingle",
        init<SphericalSingle::surface_id_type,
             SphericalSingle::particle_id_pair,
             SphericalSingle::shell_id_pair,
             SphericalSingle::reaction_rule_vector const&>())
        .add_property("particle",
            make_function(&SphericalSingle::particle,
                return_value_policy<return_by_value>()))
        .add_property("shell",
            make_function(&SphericalSingle::shell,
                return_value_policy<return_by_value>()))
        .add_property("reactions",
            make_function(&SphericalSingle::reactions,
                return_value_policy<return_by_value>()))
        .add_property("k_tot",
            make_function(&SphericalSingle::k_tot,
                return_value_policy<return_by_value>()))
        ;

    typedef _EGFRDSimulator::cylindrical_single_type CylindricalSingle;
    class_<CylindricalSingle, bases<Domain>,
           boost::shared_ptr<CylindricalSingle> >(
        "_CylindricalSingle",
        init<CylindricalSingle::surface_id_type,
             CylindricalSingle::particle_id_pair,
             CylindricalSingle::shell_id_pair,
             CylindricalSingle::reaction_rule_vector const&>())
        .add_property("particle",
            make_function(&CylindricalSingle::particle,
                return_value_policy<return_by_value>()))
        .add_property("shell",
            make_function(&CylindricalSingle::shell,
                return_value_policy<return_by_value>()))
        .add_property("reactions",
            make_function(&CylindricalSingle::reactions,
                return_value_policy<return_by_value>()))
        .add_property("k_tot",
            make_function(&CylindricalSingle::k_tot,
                return_value_policy<return_by_value>()))
        ;

    typedef _EGFRDSimulator::spherical_pair_type SphericalPair;
    class_<SphericalPair, bases<Domain>,
           boost::shared_ptr<SphericalPair> >(
        "_SphericalPair",
        init<SphericalPair::surface_id_type,
             SphericalPair::particle_id_pair,
             SphericalPair::particle_id_pair,
             SphericalPair::shell_id_pair,
             SphericalPair::length_type,
             SphericalPair::length_type>())
        .add_property("particles",
            make_function(&SphericalPair::particles,
                return_value_policy<return_by_value>()))
        .add_property("shell",
            make_function(&SphericalPair::shell,
                return_value_policy<return_by_value>()))
        .add_property("r0",
            make_function(&SphericalPair::r0,
                return_value_policy<return_by_value>()))
        .add_property("rt",
            make_function(&SphericalPair::rt,
                return_value_policy<return_by_value>()))
        .add_property("a_R",
            make_function(&SphericalPair::a_R,
                return_value_policy<return_by_value>()))
        .add_property("a_r",
            make_function(&SphericalPair::a_r,
                return_value_policy<return_by_value>()))
        .add_property("sigma",
            make_function(&SphericalPair::sigma,
                return_value_policy<return_by_value>()))
        .add_property("D_tot",
            make_function(&SphericalPair::D_tot,
                return_value_policy<return_by_value>()))
        .add_property("D_geom",
            make_function(&SphericalPair::D_geom,
                return_value_policy<return_by_value>()))
        .add_property("D_R",
            make_function(&SphericalPair::D_R,
                return_value_policy<return_by_value>()));
    peer::converters::register_range_to_tuple_converter<SphericalPair::particle_array_type>();

    typedef _EGFRDSimulator::cylindrical_pair_type CylindricalPair;
    class_<CylindricalPair, bases<Domain>,
           boost::shared_ptr<CylindricalPair> >(
        "_CylindricalPair",
        init<CylindricalPair::surface_id_type,
             CylindricalPair::particle_id_pair,
             CylindricalPair::particle_id_pair,
             CylindricalPair::shell_id_pair,
             CylindricalPair::length_type,
             CylindricalPair::length_type>())
        .add_property("particles",
            make_function(&CylindricalPair::particles,
                return_value_policy<return_by_value>()))
        .add_property("shell",
            make_function(&CylindricalPair::shell,
                return_value_policy<return_by_value>()))
        .add_property("r0",
            make_function(&CylindricalPair::r0,
                return_value_policy<return_by_value>()))
        .add_property("rt",
            make_function(&CylindricalPair::rt,
                return_value_policy<return_by_value>()))
        .add_property("a_R",
            make_function(&CylindricalPair::a_R,
                return_value_policy<return_by_value>()))
        .add_property("a_r",
            make_function(&CylindricalPair::a_r,
                return_value_policy<return_by_value>()))
        .add_property("sigma",
            make_function(&CylindricalPair::sigma,
                return_value_policy<return_by_value>()))
        .add_property("D_tot",
            make_function(&CylindricalPair::D_tot,
                return_value_policy<return_by_value>()))
        .add_property("D_geom",
            make_function(&CylindricalPair::D_geom,
                return_value_policy<return_by_value>()))
        .add_property("D_R",
            make_function(&CylindricalPair::D_R,
                return_value_policy<return_by_value>()));
    peer::converters::register_range_to_tuple_converter<CylindricalPair::particle_array_type>();

    class_<_EGFRDSimulator, boost::noncopyable>("_EGFRDSimulator",
            init<CyclicWorld&,
                 GSLRandomNumberGenerator&,
                 NetworkRulesWrapper const&>())
        .def("new_spherical_shell", &_EGFRDSimulator::new_spherical_shell)
        .def("new_cylindrical_shell", &_EGFRDSimulator::new_spherical_shell)
        .def("get_spherical_shell",
            &_EGFRDSimulator::new_spherical_shell,
            return_value_policy<return_by_value>())
        .def("get_cylindrical_shell",
            &_EGFRDSimulator::new_spherical_shell,
            return_value_policy<return_by_value>())
        .def("get_domain", &_EGFRDSimulator::get_domain)
        .def("update_domain", &_EGFRDSimulator::update_domain)
        .def("__iter__", &_EGFRDSimulator::get_domains,
                return_value_policy<return_by_value>())
        ;

    peer::wrappers::generator_wrapper<ptr_generator<_EGFRDSimulator::domain_id_pair_generator, std::auto_ptr<_EGFRDSimulator::domain_id_pair_generator> > >::__register_class("DomainIDPairGenerator");

    class_<_BDPropagator, boost::noncopyable>(
        "BDPropagator", init<
            egfrd_simulator_traits_type::world_type::particle_container_type&,
            egfrd_simulator_traits_type::network_rules_type const&,
            egfrd_simulator_traits_type::rng_type&,
            egfrd_simulator_traits_type::time_type,
            int,
            peer::wrappers::pyiterable_range<world_traits_type::particle_id_type> >())
        .def(init<
            egfrd_simulator_traits_type::world_type::particle_container_type&,
            egfrd_simulator_traits_type::network_rules_type const&,
            egfrd_simulator_traits_type::rng_type&,
            egfrd_simulator_traits_type::time_type,
            int,
            get_select_first_range<CyclicWorld::particle_id_pair_range>::type>())
        .add_property("reactions",
            peer::util::range_from_range<
                _BDPropagator::reaction_rules_range,
                _BDPropagator, &_BDPropagator::get_reactions>())
        .add_property("rejected_move_count",
            &_BDPropagator::get_rejected_move_count)
        .def("__call__", &_BDPropagator::operator())
        .def("propagate_all", &_BDPropagator_propagate_all)
        ;

    class_<MultiParticleContainer<egfrd_simulator_traits_type>,
           bases<particle_container_type>, boost::noncopyable>(
        "MultiParticleContainer", init<CyclicWorld&>());

    peer::converters::register_pyiterable_range_converter<world_traits_type::particle_id_type>();

    domain_id_pair_converter::__register();

    peer::converters::register_iterable_to_ra_container_converter<boost::array<world_traits_type::length_type, 3>, 3>();

    typedef StructureUtils<egfrd_simulator_traits_type> _StructureUtils;
    def("create_planar_surface", &_StructureUtils::create_planar_surface,
            return_value_policy<manage_new_object>());
    def("create_spherical_surface", &_StructureUtils::create_spherical_surface,
            return_value_policy<manage_new_object>());
    def("create_cylindrical_surface", &_StructureUtils::create_cylindrical_surface,
            return_value_policy<manage_new_object>());
    def("create_cuboidal_region", &_StructureUtils::create_cuboidal_region,
            return_value_policy<manage_new_object>());
    def("draw_bd_displacement", &draw_bd_displacement);

    peer::util::to_native_converter<world_traits_type::species_id_type, species_info_to_species_id_converter>();

    peer::RandomNumberGenerator<GSLRandomNumberGenerator>::__register_class("RandomNumberGenerator");
    def("create_gsl_rng", &create_gsl_rng<gsl_rng_mt19937>);
    def("create_static_gsl_rng", &create_static_gsl_rng);

    def("do_check_overlap", &do_check_overlap, return_value_policy<return_by_value>());

    peer::util::register_scalar_to_native_converters();
}
