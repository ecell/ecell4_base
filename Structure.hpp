#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <sstream>
#include "Vector3.hpp"

template<typename Ttraits_>
class Structure
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::rng_type rng_type;
    typedef typename traits_type::structure_id_type identifier_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::position_type position_type;

public:
    virtual ~Structure() {}

    identifier_type const& id() const
    {
        return id_;
    }

    virtual bool operator==(Structure const& rhs) const
    {
        return id_ == rhs.id();
    }

    bool operator!=(Structure const& rhs) const
    {
        return !operator==(rhs);
    }

    virtual position_type random_position(rng_type& rng) const = 0;

    virtual position_type random_vector(length_type const& r, rng_type& rng) const = 0;

    virtual position_type bd_displacement(length_type const& r, rng_type& rng) const = 0;

    virtual std::size_t hash() const
    {
#if defined(HAVE_TR1_FUNCTIONAL)
        using std::tr1::hash;
#elif defined(HAVE_STD_HASH)
        using std::hash;
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
        using boost::hash;
#endif
        return hash<identifier_type>()(id_);
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Structure(" << id() << ")";
        return out.str();
    }

    Structure(identifier_type const& id)
        : id_(id) {}

protected:
    identifier_type id_;
};

template<typename Tstrm, typename Ttraits, typename T_traits>
inline std::basic_ostream<Tstrm, Ttraits>& operator<<(std::basic_ostream<Tstrm, Ttraits>& strm, const Structure<T_traits>& v)
{
    strm << v.as_string(); 
    return strm;
}

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename Ttraits>
struct hash<Structure<Ttraits> >
{
    typedef Structure<Ttraits> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return val.hash();
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* STRUCTURE_HPP */
