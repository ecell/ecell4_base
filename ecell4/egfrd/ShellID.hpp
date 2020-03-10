#ifndef SHELL_ID_HPP
#define SHELL_ID_HPP

#include <functional>
#include <ostream>
// #include "Identifier.hpp"
#include <ecell4/core/Identifier.hpp>

namespace ecell4
{
namespace egfrd
{

struct ShellID: public ecell4::Identifier<ShellID, unsigned long long, int>
{
    typedef ecell4::Identifier<ShellID, unsigned long long, int> base_type;

    ShellID(value_type const& value = value_type(0, 0))
        : base_type(value) {}
};

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const ShellID& v)
{
    strm << "ShellID(" << v().first << ":" << v().second << ")";
    return strm;
}

} //egfrd
} //ecell4

namespace std {
template<>
struct hash<ecell4::egfrd::ShellID>
{
    std::size_t operator()(ecell4::egfrd::ShellID const& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};
} // std
#endif /* SHELL_ID_HPP */
