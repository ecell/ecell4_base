#ifndef ECELL_BD_MAP_UTILITUY
#define ECELL_BD_MAP_UTILITUY
#include <stdexcept>

namespace ecell4
{
namespace bd
{

template<typename mapT>
inline typename mapT::mapped_type const&
const_at(mapT const& m, typename mapT::key_type const& k)
{
    const typename mapT::const_iterator iter = m.find(k);
    if(iter == m.end()) throw std::out_of_range("ecell4::bd::const_at(mapT, key)");
    return iter->second;
}

} // bd
} // ecell4
#endif //ECELL_BD_MAP_UTILITUY
