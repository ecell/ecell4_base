#ifndef PY_HASH_SUPPORT_HPP
#define PY_HASH_SUPPORT_HPP

#include "get_mapper_mf.hpp"

template<typename Tval_>
struct get_mapper_mf<boost::python::object, Tval_>
{
#if HAVE_UNORDERED_MAP || HAVE_TR1_UNORDERED_MAP || HAVE_EXT_HASH_MAP
    struct hasher: public std::unary_function<boost::python::object, std::size_t>
    {
        typedef boost::python::object argument_type;
        typedef std::size_t result_type;

        result_type operator()(const argument_type& arg) const
        {
            return static_cast<result_type>((long)PyObject_Hash(arg.ptr()));
        }
    };
#endif


#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<boost::python::object, Tval_, hasher> type;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<boost::python::object, Tval_, hasher> type;
#elif HAVE_EXT_HASH_MAP
    typedef __gnu_cxx::hash_map<boost::python::object, Tval_, hasher> type;
#else 
    typedef std::map<boost::python::object, Tval_> type;
#endif
};

#endif /* PY_HASH_SUPPORT_HPP */
