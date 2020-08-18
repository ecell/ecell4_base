#ifndef ECELL4_NGFRD_DOMAIN_ID_HPP
#define ECELL4_NGFRD_DOMAIN_ID_HPP

#include <ecell4/core/Identifier.hpp>
#include <ostream>

namespace ecell4
{
namespace ngfrd
{

struct DomainID final : public ecell4::Identifier<DomainID, std::uint64_t, std::int32_t>
{
    using base_type = Identifier<DomainID, std::uint64_t, std::int32_t>;

    DomainID(const value_type& value = value_type(0, 0))
        : base_type(value)
    {}
};

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& strm, const DomainID& v)
{
    strm << "DomainID(" << v().first << ":" << v().second << ")";
    return strm;
}

} //ngfrd
} //ecell4

namespace std {
template<>
struct hash<ecell4::ngfrd::DomainID>
{
    std::size_t operator()(ecell4::ngfrd::DomainID const& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};
} // std
#endif// ECELL4_NGFRD_DOMAIN_ID_HPP
