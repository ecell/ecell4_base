#ifndef PEER_UTIL_EXCEPTION_TRANSLATORS_HPP
#define PEER_UTIL_EXCEPTION_TRANSLATORS_HPP

#include <Python.h>
#include <stdexcept>
#include <boost/python/exception_translator.hpp>

namespace peer { namespace util {

namespace detail
{
    template<PyObject* const& Vpytype_object, typename Texc>
    inline void exception_translator(Texc const& exc)
    {
        PyErr_SetString(Vpytype_object, exc.what());
    }
} // namespace detail

template<PyObject* const& Vpytype_object, typename Texc>
void register_exception_translator()
{
    boost::python::register_exception_translator<Texc>(
        &detail::exception_translator<Vpytype_object, Texc>);
}

inline void register_std_exception_translator()
{
    register_exception_translator<PyExc_RuntimeError, std::exception>();
    register_exception_translator<PyExc_ArithmeticError, std::domain_error>();
    register_exception_translator<PyExc_OverflowError, std::overflow_error>();
    register_exception_translator<PyExc_OverflowError, std::underflow_error>();
    register_exception_translator<PyExc_OverflowError, std::range_error>();
}

} } // namespace peer::util

#endif /* PEER_UTIL_EXCEPTION_TRANSLATORS_HPP */
