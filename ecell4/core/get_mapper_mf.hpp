#ifndef ECELL4_UTILS_GET_MAPPER_MF_HPP
#define ECELL4_UTILS_GET_MAPPER_MF_HPP

#include <ecell4/core/config.h>
#include <unordered_map>
#include <vector>

namespace ecell4
{

namespace utils
{

template<typename Tkey_, typename Tval_>
struct get_mapper_mf
{
    typedef std::unordered_map<Tkey_, Tval_> type;
};

template<typename Tmap_>
void retrieve_keys(Tmap_ map, std::vector<typename Tmap_::key_type>& keys)
{
    for (typename Tmap_::const_iterator itr(map.begin());
            itr != map.end(); ++itr)
    {
        keys.push_back(itr->first);
    }
}

} // utils

} // ecell4

#endif /* ECELL4_UTILS_GET_MAPPER_MF_HPP */
