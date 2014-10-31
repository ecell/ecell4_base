#ifndef PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP
#define PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP

#include "peer/wrappers/generator/generator_wrapper.hpp"

namespace peer { namespace converters {

template<typename Tgen_>
struct generator_to_pyiterator_converter 
{
    typedef Tgen_ native_type;
    typedef peer::wrappers::generator_wrapper<native_type> wrapper_type;

    static PyObject* convert(native_type const& impl)
    {
        return wrapper_type::create(impl);
    }

    static PyTypeObject* get_pytype()
    {
        return &wrapper_type::__class__;
    }
};

template<typename Tgen_, typename Tholder_>
struct ptr_generator_to_pyiterator_converter
{
    typedef Tgen_* native_type;
    typedef ptr_generator<Tgen_, Tholder_> ptr_generator_type;
    typedef peer::wrappers::generator_wrapper<ptr_generator_type> wrapper_type;

    static PyObject* convert(native_type const& impl)
    {
        if (impl)
        {
            Tholder_ ptr(impl);
            return wrapper_type::create(ptr_generator_type(ptr));
        }
        return boost::python::incref(Py_None);
    }

    static PyTypeObject* get_pytype()
    {
        return &wrapper_type::__class__;
    }
};

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP */
