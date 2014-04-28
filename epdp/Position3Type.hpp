#ifndef POSITION3_TRAITS_HPP
#define POSITION3_TRAITS_HPP

// This file containing the reference to the ecell4::Position3 and 
//   some template traits of ecell4::Position3.
//
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <ecell4/core/Position3.hpp>

#include <boost/array.hpp>
#include "utils/array_traits.hpp"
#include "linear_algebra.hpp"

template <std::size_t N_>
struct is_vector<ecell4::Position3, N_>: public boost::mpl::true_ {};

template <>
struct element_type_of<ecell4::Position3>
{
    typedef ecell4::Position3::value_type type;
};

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<>
struct hash<ecell4::Position3>
{
    typedef ecell4::Position3 argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<argument_type::value_type>()(val[0]) ^
            hash<argument_type::value_type>()(val[1]) ^
            hash<argument_type::value_type>()(val[2]);
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif



#endif /* POSITION3_TRAITS_HPP */
