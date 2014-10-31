#ifndef PEER_UTIL_TO_PYTHON_CONVERTER_FUN_HPP
#define PEER_UTIL_TO_PYTHON_CONVERTER_FUN_HPP

#include <functional>
#include <boost/python.hpp>

namespace peer { namespace util {

template<typename T_>
struct to_python_converter_fun
    : public std::unary_function<T_, boost::python::object>
{
    typedef T_ argument_type;
    typedef boost::python::object result_type;

    result_type operator()(argument_type const& src) const
    {
        return boost::python::object(src);
    }
};

} } // namespace peer::util

#endif /* PEER_UTIL_TO_PYTHON_CONVERTER_FUN_HPP */
