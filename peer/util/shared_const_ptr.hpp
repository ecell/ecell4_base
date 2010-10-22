#ifndef PEER_UTIL_SHARED_CONST_PTR
#define PEER_UTIL_SHARED_CONST_PTR

#include "peer/util/to_native_converter.hpp"
#include <boost/python/converter/shared_ptr_from_python.hpp>
// stolen from <boost/python/converter/shared_ptr_from_python>

namespace peer { namespace util {

namespace detail {

template<typename T_>
struct shared_const_ptr_from_python
{
    static PyTypeObject const* expected_pytype()
    {
        return boost::python::converter::expected_from_python_type_direct<T_ const>::get_pytype();
    }

    static void* convertible(PyObject* p)
    {
        if (p == Py_None)
            return p;
        
        return boost::python::converter::get_lvalue_from_python(p, boost::python::converter::registered<T_>::converters);
    }
    
    static void construct(PyObject* source, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        void* const storage = ((boost::python::converter::rvalue_from_python_storage<boost::shared_ptr<T_ const> >*)data)->storage.bytes;
        // Deal with the "None" case.
        if (data->convertible == source)
            new (storage) boost::shared_ptr<T_ const>();
        else
        {
            boost::shared_ptr<void> hold_convertible_ref_count(
              (void*)0, boost::python::converter::shared_ptr_deleter(boost::python::handle<>(boost::python::borrowed(source))) );
            // use aliasing constructor
            new (storage) boost::shared_ptr<T_ const>(
                hold_convertible_ref_count,
                static_cast<T_ const*>(data->convertible));
        }
        
        data->convertible = storage;
    }
};

} // namespace detail

template<typename T>
void register_shared_const_ptr_from_python()
{
    peer::util::to_native_converter<boost::shared_ptr<T const>,
        detail::shared_const_ptr_from_python<T> >();
}

} }

#endif /* PEER_UTIL_SHARED_CONST_PTR */
