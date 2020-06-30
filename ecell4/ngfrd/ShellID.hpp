#ifndef ECELL4_NGFRD_SHELL_ID_HPP
#define ECELL4_NGFRD_SHELL_ID_HPP

#include <ostream>
#include <ecell4/core/Identifier.hpp>

namespace ecell4
{
namespace ngfrd
{

struct ShellID final : public ecell4::Identifier<ShellID, std::uint64_t, std::int32_t>
{
    using base_type = Identifier<ShellID, std::uint64_t, std::int32_t>;

    ShellID(const value_type& value = value_type(0, 0))
        : base_type(value)
    {}
};

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& strm, const ShellID& v)
{
    strm << "ShellID(" << v().first << ":" << v().second << ")";
    return strm;
}

} //ngfrd
} //ecell4

namespace std {
template<>
struct hash<ecell4::ngfrd::ShellID>
{
    std::size_t operator()(ecell4::ngfrd::ShellID const& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};
} // std
#endif /* SHELL_ID_HPP */
