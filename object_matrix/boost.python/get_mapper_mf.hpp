#ifndef GET_MAPPER_MF_HPP
#define GET_MAPPER_MF_HPP

#if HAVE_UNORDERED_MAP
#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#elif HAVE_EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <map>
#endif /* HAVE_UNORDERED_MAP */

template<typename Tkey_, typename Tval_>
struct get_mapper_mf
{
#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<Tkey_, Tval_> type;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<Tkey_, Tval_> type;
#elif HAVE_EXT_HASH_MAP
    typedef __gnu_cxx::hash_map<Tkey_, Tval_> type;
#else 
    typedef std::map<Tkey_, Tval_> type;
#endif
};

#endif /* GET_MAPPER_MF_HPP */
