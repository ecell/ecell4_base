#ifndef ECELL4_EGFRD_DOMAIN_ID_HPP
#define ECELL4_EGFRD_DOMAIN_ID_HPP

#include <ostream>
#include <functional>
// #include "Identifier.hpp"
#include <ecell4/core/Identifier.hpp>

namespace ecell4
{
namespace egfrd
{

struct DomainID: public ecell4::Identifier<DomainID, unsigned long long, int>
{
    typedef ecell4::Identifier<DomainID, unsigned long long, int> base_type;

    DomainID(value_type const& value = value_type(0, 0))
        : base_type(value) {}
};

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const DomainID& v)
{
    strm << "DID(" << v().first << ":" << v().second << ")";
    return strm;
}

} // egfrd
} // ecell4

namespace std {

template<>
struct hash<ecell4::egfrd::DomainID>
{
    std::size_t operator()(ecell4::egfrd::DomainID const& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};

} // std
#endif /* DOMAIN_ID_HPP */
