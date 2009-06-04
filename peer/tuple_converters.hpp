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
                Ttuple_, detail::tuple_to_tuple_converter<Ttuple_> >();
            boost::python::to_python_converter<
                boost::shared_ptr<Ttuple_>,
                detail::tuple_to_tuple_converter<boost::shared_ptr<Ttuple_> > >();
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
