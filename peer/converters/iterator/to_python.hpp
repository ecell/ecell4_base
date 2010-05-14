#ifndef PEER_CONVERTERS_ITERATOR_TO_PYTHON_HPP
#define PEER_CONVERTERS_ITERATOR_TO_PYTHON_HPP

#include <boost/python.hpp>
#include "peer/wrappers/iterator/stl_iterator_wrapper.hpp"

namespace peer { namespace converters {

template<typename Trange_, typename Tholder_ = void*, typename Trcg_ = boost::python::return_by_value>
struct stl_iterator_range_converter
{
    typedef Trange_ native_type;

    static PyObject* convert(native_type const& v)
    {
        return reinterpret_cast<PyObject*>(stl_iterator_wrapper<typename boost::range_const_iterator<native_type>::type, Tholder_, Trcg_>::create(v));
    }
};

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_ITERATOR_TO_PYTHON_HPP */
