#ifndef PEER_CONVERTERS_TUPLE_HPP
#define PEER_CONVERTERS_TUPLE_HPP

#include <boost/shared_ptr.hpp>
#include "peer/converters/tuple/to_python.hpp"
#include "peer/converters/tuple/from_python.hpp"
#include "peer/utils.hpp"

namespace peer { namespace converters {

template<typename Ttuple_>
inline void register_tuple_converter()
{
    static bool registered = false;
    if (!registered)
    {
        boost::python::to_python_converter<
            Ttuple_, tuple_to_pytuple_converter<Ttuple_> >();
        boost::python::to_python_converter<
            boost::shared_ptr<Ttuple_>,
            tuple_to_pytuple_converter<boost::shared_ptr<Ttuple_> > >();
        peer::util::to_native_converter<Ttuple_, pytuple_to_tuple_converter<Ttuple_> >();
        registered = true;
    }
}

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_TUPLE_HPP */
