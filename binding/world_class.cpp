#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "World.hpp"
#include "binding_common.hpp"
#include "peer/converters/sequence.hpp"

namespace binding {

template<typename Timpl_>
struct particle_id_pair_range_select_first_converter
{
    typedef Timpl_ native_type;

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


void register_world_class()
{
    register_world_class<World, ParticleContainer, EGFRDSimulator>("World");
    particle_id_pair_range_select_first_converter<
        get_select_first_range<
            World::particle_id_pair_range>::type>::__register();
    peer::converters::register_range_to_tuple_converter<twofold_container<ParticleID> >();
    peer::converters::register_iterable_to_range_converter<twofold_container<ParticleID> >();
}

} // namesapce binding
