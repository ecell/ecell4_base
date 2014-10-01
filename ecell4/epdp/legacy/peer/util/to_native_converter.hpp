#ifndef PEER_UTIL_TO_NATIVE_CONVERTER_HPP
#define PEER_UTIL_TO_NATIVE_CONVERTER_HPP

#include <boost/python.hpp>

namespace peer { namespace util {

template<typename Tnative_, typename Tconverter_>
inline void to_native_converter()
{
    boost::python::converter::registry::push_back(
            &Tconverter_::convertible,
            reinterpret_cast<
                    boost::python::converter::constructor_function>(
                        &Tconverter_::construct),
            boost::python::type_id<Tnative_>());
}

template<typename Tnative_, typename Tconverter_>
inline void to_native_lvalue_converter()
{
    boost::python::converter::registry::insert(
            &Tconverter_::convert,
            boost::python::type_id<Tnative_>(),
            &Tconverter_::expected_pytype);
}

} } // namespace peer::util

#endif /* PEER_UTIL_TO_NATIVE_CONVERTER_HPP */
