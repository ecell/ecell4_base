#ifndef GET_MAPPER_MF_HPP
#define GET_MAPPER_MF_HPP

#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(HAVE_BOOST_UNORDERED_MAP_HPP)
#include <boost/unordered_map.hpp>
#else
#include <map>
#endif /* HAVE_UNORDERED_MAP */

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

#endif /* GET_MAPPER_MF_HPP */
