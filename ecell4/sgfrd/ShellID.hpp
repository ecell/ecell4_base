#ifndef ECELL4_SGFRD_SHELL_ID
#define ECELL4_SGFRD_SHELL_ID
#include <ecell4/core/config.h>
#include <ecell4/core/hash.hpp>
#include <ecell4/core/Identifier.hpp>

#include <ostream>

struct ShellID: public ecell4::Identifier<ShellID, unsigned long long, int>
{
    typedef ecell4::Identifier<ShellID, unsigned long long, int> base_type;

    ShellID(value_type const& value = value_type(0, 0))
        : base_type(value) {}
};

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<ShellID>
{
    std::size_t operator()(ShellID const& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};

ECELL4_DEFINE_HASH_END()

template<typename charT, typename traitsT>
inline std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& strm, const ShellID& v)
{
    strm << "ShellID(" << v().first << ":" << v().second << ")";
    return strm;
}

#endif /* ECELL4_SGFRD_SHELL_ID */
