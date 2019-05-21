#ifndef ECELL4_HASH_HPP
#define ECELL4_HASH_HPP

#if defined(HAVE_STD_HASH)

#include <functional>
#define ECELL4_DEFINE_HASH_BEGIN() namespace std {
#define ECELL4_DEFINE_HASH_END() } // std
#define ECELL4_HASH_STRUCT std::hash

#elif defined(HAVE_TR1_FUNCTIONAL)

#include <tr1/functional>
#define ECELL4_DEFINE_HASH_BEGIN() namespace std { namespace tr1 {
#define ECELL4_DEFINE_HASH_END() } } // tr1 std
#define ECELL4_HASH_STRUCT std::tr1::hash

#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)

#include <boost/functional/hash.hpp>
#define ECELL4_DEFINE_HASH_BEGIN() namespace boost {
#define ECELL4_DEFINE_HASH_END() } // boost
#define ECELL4_HASH_STRUCT boost::hash

#else

#define ECELL4_DEFINE_HASH_BEGIN()
#define ECELL4_DEFINE_HASH_END()
#define ECELL4_HASH_STRUCT hash

#endif

#endif /* ECELL4_HASH_HPP */
