#ifndef ECELL4_SGFRD_DOMAIN_ID
#define ECELL4_SGFRD_DOMAIN_ID
#include <ecell4/core/config.h>
#include <ecell4/core/hash.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ostream>

namespace ecell4
{
namespace sgfrd
{
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

template<typename charT, typename traitsT>
inline std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& strm, const DomainID& v)
{
    strm << "DomainID(" << v().first << ":" << v().second << ")";
    return strm;
}

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_DOMAIN_ID */
