#ifndef ECELL4_UTILS_GET_MAPPER_MF_HPP
#define ECELL4_UTILS_GET_MAPPER_MF_HPP

#include <ecell4/core/config.h>

#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(HAVE_BOOST_UNORDERED_MAP_HPP)
#include <boost/unordered_map.hpp>
#else
#include <map>
#endif /* HAVE_UNORDERED_MAP */

#include <vector>


namespace ecell4
{

namespace utils
{

/**
   a metafunction for generating a type for an efficient map algorithm.
   in the current version of C++, the following line is not accepted:
   template<typename Tkey_, typename Tval_>
   typdef std::map<Tkey_, Tval_> map_type;
   see http://www.boost.org/community/generic_programming.html#type_generator
 */
template<typename Tkey_, typename Tval_>
struct get_mapper_mf
{
#if defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_map<Tkey_, Tval_> type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_map<Tkey_, Tval_> type;
#elif defined(HAVE_BOOST_UNORDERED_MAP_HPP)
    typedef boost::unordered_map<Tkey_, Tval_> type;
#else
    typedef std::map<Tkey_, Tval_> type;
#endif
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
