#ifndef PEER_RANDOMNUMBERGENERATOR_HPP
#define PEER_RANDOMNUMBERGENERATOR_HPP

#include <boost/python.hpp>
#include "peer/compat.h"

namespace peer {

template<typename Timpl_>
struct RandomNumberGenerator
{
    typedef Timpl_ impl_type;

    static PyObject* normal(impl_type& impl, Real loc, Real scale, boost::python::object size)
    {
        if (size.ptr() == Py_None) {
            return boost::python::incref(boost::python::object(
                impl.normal(loc, scale)).ptr());
        } else {
            ssize_t len = 0;
            std::size_t num_samples = 1;
            PyObject* retval = 0;
            len = PySequence_Size(size.ptr());
            // cancel the error that would occur if size wasn't be a sequence
            if (PyErr_Occurred())
            {
                PyErr_Clear();
            }
            if (len == -1 && PyNumber_Check(size.ptr()))
            {
                npy_intp dims[1] = { PyNumber_AsSsize_t(size.ptr(), NULL) };
                if (dims[0] < 0)
                {
                    PyErr_SetString(PyExc_ValueError, "expected size >= 0");
                    boost::python::throw_error_already_set();
                }
                num_samples = dims[0];
                retval = PyArray_New(
                        &PyArray_Type, 1, dims,
                        peer::util::get_numpy_typecode<Real>::value,
                        NULL, NULL, 0, NPY_CARRAY, NULL);
            }
            else if (len >= 0)
            {
                if (len > NPY_MAXDIMS)
                {
                    PyErr_Format(PyExc_ValueError, "sequence too large; must be smaller than %d", NPY_MAXDIMS);
                    boost::python::throw_error_already_set();
                }
                npy_intp dims[len];
                for (ssize_t i = 0; i < len; ++i)
                {
                    dims[i] = boost::python::extract<npy_intp>(size[i]);
                    num_samples *= dims[i];
                }
                retval = PyArray_New(
                        &PyArray_Type, len, dims,
                        peer::util::get_numpy_typecode<Real>::value,
                        NULL, NULL, 0, NPY_CARRAY, NULL);
            }
            Real* data(reinterpret_cast<Real*>(PyArray_DATA(retval)));
            for (std::size_t i = 0; i < num_samples; ++i)
            {
                data[i] = impl.normal(loc, scale);
            }
            return retval;
        }
    }

    static void __register_class(const char* name)
    {
        using namespace boost::python;
        class_<impl_type>(name, no_init)
            .def("normal", &RandomNumberGenerator::normal)
            .def("normal", &impl_type::normal)
            .def("seed", &impl_type::seed)
            .def("get_raw", &impl_type::get_raw)
            .def("uniform", &impl_type::uniform)
            .def("uniform_int", &impl_type::uniform_int)
            .def("__call__", &impl_type::operator())
            ;
    }
};

} // namespace peer

#endif /* PEER_RANDOMNUMBERGENERATOR_HPP */
