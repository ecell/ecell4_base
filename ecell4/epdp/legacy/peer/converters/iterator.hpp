#ifndef PEER_CONVERTERS_ITERATOR_HPP
#define PEER_CONVERTERS_ITERATOR_HPP

#include "peer/converters/iterator/to_python.hpp"

namespace peer { namespace converters {

template<typename Trange, typename Tholder, typename Trcg>
inline void register_stl_iterator_range_converter()
{
    boost::python::to_python_converter<Trange, stl_iterator_range_converter<Trange, Tholder, Trcg> >();
}

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_ITERATOR_HPP */
