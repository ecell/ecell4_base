#ifndef ASSOC_CONTAINER_TRAITS_HPP
#define ASSOC_CONTAINER_TRAITS_HPP

#include <boost/range/value_type.hpp>

template<typename Tassoc_>
struct assoc_key
{
    typedef typename Tassoc_::key_type type;
};

template<typename Tassoc_>
struct assoc_value
{
    typedef typename boost::range_value<Tassoc_>::type type;
};
 
template<typename Tassoc_>
struct assoc_mapped
{
    typedef typename Tassoc_::mapped_type type;
};

#endif /* ASSOC_CONTAINER_TRAITS_HPP */
