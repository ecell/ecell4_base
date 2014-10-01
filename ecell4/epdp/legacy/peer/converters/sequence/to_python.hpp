#ifndef PEER_CONVERTERS_SEQUENCE_TO_PYTHON_HPP
#define PEER_CONVERTERS_SEQUENCE_TO_PYTHON_HPP

#include <Python.h>
#include <tupleobject.h>
#include <listobject.h>
#include <boost/python.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

namespace peer { namespace converters {

template<typename Trange_, typename Tpolicy_>
struct range_to_pyseq_converter
{
    typedef Trange_ native_type;
    typedef Tpolicy_ policy_type;
    
    static PyObject* convert(const native_type& p)
    {
        using namespace boost::python;
        PyObject* retval = policy_type::create(boost::size(p));
        Py_ssize_t idx = 0;
        for (typename boost::range_const_iterator<native_type>::type i(boost::begin(p)), e(boost::end(p)); i != e; ++i, ++idx)
        {
            policy_type::set(retval, idx, incref(object(*i).ptr()));
        }
        return retval;
    }
};

struct tuple_policy
{
    static PyObject* create(Py_ssize_t size)
    {
        return PyTuple_New(size);
    }

    static void set(PyObject* tuple, Py_ssize_t idx, PyObject* obj)
    {
        PyTuple_SET_ITEM(tuple, idx, obj);
    }
};

struct list_policy
{
    static PyObject* create(Py_ssize_t size)
    {
        return PyList_New(size);
    }

    static void set(PyObject* tuple, Py_ssize_t idx, PyObject* obj)
    {
        PyList_SET_ITEM(tuple, idx, obj);
    }
};

template<typename Trange_>
struct range_to_pytuple_converter
    : public range_to_pyseq_converter<Trange_, tuple_policy>
{
    typedef Trange_ native_type;
};

template<typename Trange_>
struct range_to_pylist_converter
    : public range_to_pyseq_converter<Trange_, list_policy>
{
    typedef Trange_ native_type;
};

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_SEQUENCE_TO_PYTHON_HPP */
