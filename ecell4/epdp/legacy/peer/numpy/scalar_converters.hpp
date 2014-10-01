#ifndef PEER_NUMPY_SCALAR_CONVERTERS_HPP
#define PEER_NUMPY_SCALAR_CONVERTERS_HPP

#include <stdexcept>
#include <complex>
#include <vector>
#include <boost/python.hpp>
#include <boost/multi_array.hpp>
#include <boost/format.hpp>
#include <numpy/arrayobject.h>
#include "peer/utils.hpp"
#include "peer/numpy/type_mappings.hpp"

namespace peer {

namespace util
{
    namespace detail
    {
        template<typename T_>
        struct scalar_to_native_converter
        {
            typedef T_ native_type;

            static void* convertible(PyObject* pyo)
            {
                return PyArray_CheckScalar(pyo) ? pyo: 0;
            }

            static void construct(PyObject* pyo,
                                  boost::python::converter::rvalue_from_python_stage1_data* data)
            {
                void* storage(reinterpret_cast<
                    boost::python::converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
                PyArray_Descr* descr = PyArray_DescrFromType(get_numpy_typecode<T_>::value);
                if (PyArray_CastScalarToCtype(pyo, storage, descr))
                {
                    PyErr_SetString(PyExc_TypeError,
                            (boost::format("Failed to cast %s to %s") % pyo->ob_type->tp_name % descr->typeobj->tp_name).str().c_str());
                    Py_DECREF(reinterpret_cast<PyObject*>(descr));
                    boost::python::throw_error_already_set();
                }
                Py_DECREF(reinterpret_cast<PyObject*>(descr));
                data->convertible = storage;
            }
        };
    } // namespace detail

    template<typename T_>
    inline void register_scalar_to_native_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            to_native_converter<T_, detail::scalar_to_native_converter<T_> >();
            registered = true;
        }
    }

    inline void register_scalar_to_native_converters()
    {
        register_scalar_to_native_converter<bool>();
        register_scalar_to_native_converter<npy_byte>();
        register_scalar_to_native_converter<npy_ubyte>();
        register_scalar_to_native_converter<npy_short>();
        register_scalar_to_native_converter<npy_ushort>();
        register_scalar_to_native_converter<npy_int>();
        register_scalar_to_native_converter<npy_uint>();
        register_scalar_to_native_converter<npy_long>();
        register_scalar_to_native_converter<npy_ulong>();
        register_scalar_to_native_converter<npy_longlong>();
        register_scalar_to_native_converter<npy_ulonglong>();
        register_scalar_to_native_converter<npy_float>();
        register_scalar_to_native_converter<npy_double>();
        register_scalar_to_native_converter<npy_longdouble>();
        register_scalar_to_native_converter<char>();
    }
} // namespace util

} // namespace peer

#endif /* PPER_NUMPY_SCALAR_CONVERTERS_HPP */

