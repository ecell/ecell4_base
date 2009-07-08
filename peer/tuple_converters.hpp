#ifndef OBJECTMATRIX_PEER_TUPLE_CONVERTERS_HPP
#define OBJECTMATRIX_PEER_TUPLE_CONVERTERS_HPP

#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/range/size.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>

#include "utils.hpp"

namespace peer {

namespace util
{
    namespace detail
    {
        template<typename Ttcell_>
        inline void build_pytuple_from_tuple(PyObject* pyt, const Ttcell_& cell,
                Py_ssize_t idx = 0)
        {
            PyTuple_SetItem(pyt, idx,
                boost::python::incref(
                    boost::python::object(cell.get_head()).ptr()));
            build_pytuple_from_tuple(pyt, cell.get_tail(), idx + 1);
        }

        template<>
        inline void build_pytuple_from_tuple<boost::tuples::null_type>(
                PyObject*, const boost::tuples::null_type&, Py_ssize_t) {}

        template<typename Ttuple_>
        struct tuple_to_pytuple_converter
        {
            typedef Ttuple_ argument_value_type;
            typedef const argument_value_type& argument_type;
            static PyObject* convert(argument_type val)
            {
                PyObject* retval =
                    PyTuple_New(boost::tuples::length<Ttuple_>::value);
                build_pytuple_from_tuple(retval, val);
                return retval;
            }
        };

        template<typename Tfirst_, typename Tsecond_>
        struct tuple_to_pytuple_converter<std::pair<Tfirst_, Tsecond_> >
        {
            typedef std::pair<Tfirst_, Tsecond_> argument_value_type;
            typedef const argument_value_type& argument_type;

            static PyObject* convert(argument_type val)
            {
                return boost::python::incref(
                        boost::python::make_tuple(
                                val.first, val.second).ptr());
            }
        };

        template<typename Ttuple_>
        struct tuple_to_pytuple_converter<boost::shared_ptr<Ttuple_> >
        {
            typedef Ttuple_ argument_value_type;
            typedef boost::shared_ptr<argument_value_type> argument_type;
            static PyObject* convert(argument_type val)
            {
                return tuple_to_pytuple_converter<argument_value_type>::convert(
                        *val);
            }
        };

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
        boost::python::extract<BOOST_PP_CAT(__v__,__n__)>(PySequence_GetItem(pyo, __n__))

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
                    boost::python::extract<Tfirst_>(PySequence_GetItem(pyo, 0)),
                    boost::python::extract<Tsecond_>(PySequence_GetItem(pyo, 1)));
            }
        };

        template<typename Trange_>
        struct range_to_pytuple_converter
        {
            typedef Trange_ native_type;
            
            static PyObject* convert(const native_type& p)
            {
                using namespace boost::python;
                PyObject* retval = PyTuple_New(boost::size(p));
                Py_ssize_t idx = 0;
                for (typename boost::range_const_iterator<native_type>::type i(boost::begin(p)), e(boost::end(p)); i != e; ++i, ++idx)
                {
                    PyTuple_SetItem(retval, idx, incref(object(*i).ptr()));
                }
                return retval;
            }
        };
    } // namespace detail

    template<typename Ttuple_>
    void register_tuple_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            boost::python::to_python_converter<
                Ttuple_, detail::tuple_to_pytuple_converter<Ttuple_> >();
            boost::python::to_python_converter<
                boost::shared_ptr<Ttuple_>,
                detail::tuple_to_pytuple_converter<boost::shared_ptr<Ttuple_> > >();
            to_native_converter<Ttuple_, detail::pytuple_to_tuple_converter<Ttuple_> >();
            registered = true;
        }
    }


    template<typename Trange_>
    void register_range_to_tuple_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            boost::python::to_python_converter<
                Trange_, detail::range_to_pytuple_converter<Trange_> >();
            registered = true;
        }
    }

} // namespace util

} // namespace peer

#endif /* OBJECTMATRIX_PEER_TUPLE_CONVERTERS_HPP */
