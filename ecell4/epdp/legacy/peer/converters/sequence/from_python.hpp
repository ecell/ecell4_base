#ifndef PEER_CONVERTERS_SEQUENCE_FROM_PYTHON_HPP
#define PEER_CONVERTERS_SEQUENCE_FROM_PYTHON_HPP

#include <boost/python.hpp>
#include <boost/range/value_type.hpp>
#include "peer/wrappers/range/pyiterable_range.hpp"

namespace peer { namespace converters {

template<typename Tcntnr_>
struct pyiterable_to_container_converter
{
    typedef Tcntnr_ native_type;

    static void* convertible(PyObject* pyo)
    {
        PyObject* const retval(PyObject_GetIter(pyo));
        if (!retval)
        {
            PyErr_Clear();
            return 0;
        }
        return retval;
    }

    static void construct(PyObject* pyo,
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage(reinterpret_cast<
            boost::python::converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
        boost::python::handle<> iter(
                reinterpret_cast<PyObject*>(data->convertible));

        data->convertible = new (storage) native_type();
        native_type& retval(*reinterpret_cast<native_type*>(data->convertible));
        for (;;)
        {
            boost::python::handle<> i(
                    boost::python::allow_null(PyIter_Next(iter.get())));
            if (!i)
            {
                if (PyErr_Occurred())
                {
                    boost::python::throw_error_already_set();
                }
                break;
            }
            retval.insert(boost::end(retval),
                boost::python::extract<
                    typename boost::range_value<native_type>::type>(
                        i.get())());
        }
    }
};

template<typename Trange_, std::size_t N_>
struct pyiterable_to_ra_container_converter
{
    typedef Trange_ native_type;

    static void* convertible(PyObject* pyo)
    {
        PyObject* const retval(PyObject_GetIter(pyo));
        if (!retval)
        {
            PyErr_Clear();
            return 0;
        }
        return retval;
    }

    static void construct(PyObject* pyo,
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage(reinterpret_cast<
            boost::python::converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
        boost::python::handle<> iter(
                reinterpret_cast<PyObject*>(data->convertible));

        data->convertible = new (storage) native_type();
        native_type& retval(*reinterpret_cast<native_type*>(data->convertible));
        std::size_t idx(0);
        for (;;)
        {
            boost::python::handle<> i(boost::python::allow_null(PyIter_Next(iter.get())));
            if (!i)
            {
                if (PyErr_Occurred())
                {
                    boost::python::throw_error_already_set();
                }
                break;
            }
            if (idx >= N_)
            {
                PyErr_Format(PyExc_ValueError, "iterable generated more than %zd items", N_);
                boost::python::throw_error_already_set();
            }
            retval[idx++] = boost::python::extract<
                    typename boost::range_value<native_type>::type>(
                        i.get())();
        }
        if (idx < N_)
        {
            PyErr_Format(PyExc_ValueError, "iterable generated less than %zd items", N_);
            boost::python::throw_error_already_set();
        }
    }
};

template<typename Tvalue_>
struct pyiterable_range_converter
{
    typedef peer::wrappers::pyiterable_range<Tvalue_> native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyType_HasFeature(Py_TYPE(pyo), Py_TPFLAGS_HAVE_ITER))
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
            converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
        data->convertible = new (storage) native_type(
            boost::python::object(boost::python::borrowed(pyo)));
    }
};


} } // namespace peer::converters

#endif /* PEER_CONVERTERS_SEQUENCE_FROM_PYTHON_HPP */
