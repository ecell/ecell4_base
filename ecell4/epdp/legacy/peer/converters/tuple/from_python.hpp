#ifndef PEER_CONVERTERS_TUPLE_FROM_PYTHON_HPP
#define PEER_CONVERTERS_TUPLE_FROM_PYTHON_HPP

#include <utility>
#include <boost/tuple/tuple.hpp>
#include <boost/python/tuple.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>

namespace peer { namespace converters {

template<typename Ttuple_>
struct pytuple_to_tuple_converter
{
    static void* convertible(PyObject* pyo)
    {
        return 0;
    }

    static void construct(PyObject* pyo,
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
    }
};

#define TUPLE_CONVERTERS_PYTUPLE_TO_TUPLE_TEMPLATE_EXTRACT(__z__, __n__, __v__) \
(boost::python::extract<BOOST_PP_CAT(__v__,__n__)>(boost::python::handle<>(PySequence_GetItem(pyo, __n__)).get()))()

#define TUPLE_CONVERTERS_PYTUPLE_TO_TUPLE_TEMPLATE(__z__, __n__, __v__) \
template<BOOST_PP_ENUM_PARAMS(__n__, typename T)> \
struct pytuple_to_tuple_converter<boost::tuples::tuple<BOOST_PP_ENUM_PARAMS(__n__, T)> > \
{ \
    typedef boost::tuples::tuple<BOOST_PP_ENUM_PARAMS(__n__, T)> native_type; \
\
    static void* convertible(PyObject* pyo) \
    { \
        if (!PySequence_Check(pyo)) \
        { \
            return 0; \
        } \
\
        if (PySequence_Size(pyo) != __n__) \
        { \
            return 0; \
        } \
        return pyo; \
    } \
\
    static void construct(PyObject* pyo, \
                          boost::python::converter::rvalue_from_python_stage1_data* data) \
    { \
        void* storage(reinterpret_cast< \
            boost::python::converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes); \
        data->convertible = new (storage) native_type( \
            BOOST_PP_ENUM(__n__, TUPLE_CONVERTERS_PYTUPLE_TO_TUPLE_TEMPLATE_EXTRACT, T)); \
    } \
};

BOOST_PP_REPEAT_FROM_TO(1, 10, TUPLE_CONVERTERS_PYTUPLE_TO_TUPLE_TEMPLATE, );

#undef TUPLE_CONVERTERS_PYTUPLE_TO_TUPLE_TEMPLATE_NULLTYPE
#undef TUPLE_CONVERTERS_PYTUPLE_TO_TUPLE_TEMPLATE_EXTRACT
#undef TUPLE_CONVERTERS_PYTUPLE_TO_TUPLE_TEMPLATE

template<typename Tfirst_, typename Tsecond_>
struct pytuple_to_tuple_converter<std::pair<Tfirst_, Tsecond_> >
{
    typedef std::pair<Tfirst_, Tsecond_> native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PySequence_Check(pyo))
        {
            return 0;
        }

        if (PySequence_Size(pyo) != 2)
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo,
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage(reinterpret_cast<
            boost::python::converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
        data->convertible = new (storage) native_type(
            boost::python::extract<Tfirst_>(boost::python::handle<>(PySequence_GetItem(pyo, 0)).get())(),
            boost::python::extract<Tsecond_>(boost::python::handle<>(PySequence_GetItem(pyo, 1)).get())());
    }
};

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_TUPLE_FROM_PYTOHN_HPP */
