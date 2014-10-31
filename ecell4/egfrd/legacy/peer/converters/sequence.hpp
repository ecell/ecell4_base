#ifndef PEER_CONVERTERS_SEQUENCE_HPP
#define PEER_CONVERTERS_SEQUENCE_HPP

#include "peer/wrappers/range/pyiterable_range.hpp"
#include "peer/converters/sequence/from_python.hpp"
#include "peer/converters/sequence/to_python.hpp"
#include "peer/utils.hpp"

namespace peer { namespace converters {

template<typename Trange_>
inline void register_range_to_tuple_converter()
{
    static bool registered = false;
    if (!registered)
    {
        boost::python::to_python_converter<
            Trange_, range_to_pytuple_converter<Trange_> >();
        registered = true;
    }
}

template<typename Trange_>
inline void register_range_to_list_converter()
{
    static bool registered = false;
    if (!registered)
    {
        boost::python::to_python_converter<
            Trange_, range_to_pylist_converter<Trange_> >();
        registered = true;
    }
}

template<typename Trange_>
inline void register_iterable_to_range_converter()
{
    static bool registered = false;
    if (!registered)
    {
        peer::util::to_native_converter<Trange_, pyiterable_to_container_converter<Trange_> >();
        registered = true;
    }
}

template<typename Trange_, std::size_t N_>
inline void register_iterable_to_ra_container_converter()
{
    static bool registered = false;
    if (!registered)
    {
        peer::util::to_native_converter<Trange_, pyiterable_to_ra_container_converter<Trange_, N_> >();
        registered = true;
    }
}

template<typename Tvalue_>
inline void register_pyiterable_range_converter()
{
    static bool registered = false;
    if (!registered)
    {
        peer::util::to_native_converter<
            peer::wrappers::pyiterable_range<Tvalue_>,
            peer::converters::pyiterable_range_converter<Tvalue_> >();
        registered = true;
    }
}

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_SEQUENCE_HPP */
