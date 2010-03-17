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

template<typename Tid_>
class Structure
{
public:
    typedef Tid_ identifier_type;

public:
    virtual ~Structure() {}

    identifier_type const& id() const
    {
        return id_;
    }

    bool operator==(Structure const& rhs) const
    {
        return id_ == rhs.id();
    }

    bool operator!=(Structure const& rhs) const
    {
        return !operator==(rhs);
    }

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

template<typename Tstrm_, typename Ttraits_, typename Tid_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Structure<Tid_>& v)
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

template<typename Tid_>
struct hash<Structure<Tid_> >
{
    typedef Structure<Tid_> argument_type;

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
