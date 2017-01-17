#ifndef ECELL_BD_MAP_UTILITUY
#define ECELL_BD_MAP_UTILITUY

namespace ecell4
{
namespace bd
{

// its nice to have static_assert or enable_if
template<typename mapT>
inline typename mapT::mapped_type const&
const_at(mapT const& m, typename mapT::key_type const& k)
{
    return m.find(k)->second;
}

} // bd
} // ecell4
#endif //ECELL_BD_MAP_UTILITUY
