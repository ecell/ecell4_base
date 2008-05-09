#ifndef OBJECTMATRIX_PEER_TUPLE_CONVERTERS_HPP
#define OBJECTMATRIX_PEER_TUPLE_CONVERTERS_HPP

#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

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
        struct tuple_to_tuple_converter
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
        struct tuple_to_tuple_converter<std::pair<Tfirst_, Tsecond_> >
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
        struct tuple_to_tuple_converter<boost::shared_ptr<Ttuple_> >
        {
            typedef Ttuple_ argument_value_type;
            typedef boost::shared_ptr<argument_value_type> argument_type;
            static PyObject* convert(argument_type val)
            {
                return tuple_to_tuple_converter<argument_value_type>::convert(
                        *val);
            }
        };
    } // namespace detail

    template<typename Ttuple_>
    void register_tuple_converter()
    {
        boost::python::to_python_converter<
            Ttuple_, detail::tuple_to_tuple_converter<Ttuple_> >();
        boost::python::to_python_converter<
            boost::shared_ptr<Ttuple_>,
            detail::tuple_to_tuple_converter<boost::shared_ptr<Ttuple_> > >();
    }

} // namespace util

} // namespace peer

#endif /* OBJECTMATRIX_PEER_TUPLE_CONVERTERS_HPP */
