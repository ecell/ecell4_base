#ifndef SPECIES_TYPE_ID_HPP
#define SPECIES_TYPE_ID_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif HAVE_STD_HASH
#include <functional>
#elif defined(BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif
#include "identifier.hpp"

struct species_type_id: public identifier<species_type_id, unsigned long long, int>
{
    typedef identifier<species_type_id, unsigned long long, int> base_type;

    species_type_id(value_type const& value = value_type(0, 0))
        : base_type(value) {}
};

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif HAVE_STD_HASH
namespace std {
#elif defined(BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<>
struct hash<species_type_id>
{
    std::size_t operator()(species_type_id const& val)
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif HAVE_STD_HASH
} // namespace std
#elif defined(BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const species_type_id& v)
{
    strm << "species_type_id(" << v().first << ":" << v().second << ")";
    return strm;
}

#endif /* SPECIES_TYPE_ID_HPP */
