#ifndef GET_DEFAULT_IMPL_HPP
#define GET_DEFAULT_IMPL_HPP

#include <map>
#include <set>
#include <vector>

#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(HAVE_BOOST_UNORDERED_MAP_HPP)
#include <boost/unordered_map.hpp>
#endif /* HAVE_UNORDERED_MAP */

namespace get_default_impl
{
    namespace std
    {
        template<typename Tkey_, typename Tval_>
        struct map
        {
            typedef ::std::map<Tkey_, Tval_> type;
        };

        template<typename Tkey_, typename Tval_>
        struct unordered_map
        {
#if defined(HAVE_UNORDERED_MAP)
            typedef ::std::unordered_map<Tkey_, Tval_> type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
            typedef ::std::tr1::unordered_map<Tkey_, Tval_> type;
#elif defined(HAVE_BOOST_UNORDERED_MAP_HPP)
            typedef ::boost::unordered_map<Tkey_, Tval_> type;
#endif
        };

        template<typename Tval_>
        struct vector
        {
            typedef ::std::vector<Tval_> type; 
        };
    } // std
} // namespace get_default_impl

#endif /* GET_DEFAULT_IMPL_HPP */
