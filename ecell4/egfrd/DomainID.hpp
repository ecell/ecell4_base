#ifndef DOMAIN_ID_HPP
#define DOMAIN_ID_HPP

#include <ecell4/core/hash.hpp>

#include <ostream>
// #include "Identifier.hpp"
#include <ecell4/core/Identifier.hpp>

struct DomainID: public ecell4::Identifier<DomainID, unsigned long long, int>
{
    typedef ecell4::Identifier<DomainID, unsigned long long, int> base_type;

    DomainID(value_type const& value = value_type(0, 0))
        : base_type(value) {}
};

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<DomainID>
{
    std::size_t operator()(DomainID const& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};

ECELL4_DEFINE_HASH_END()

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const DomainID& v)
{
    strm << "DID(" << v().first << ":" << v().second << ")";
    return strm;
}

#endif /* DOMAIN_ID_HPP */
