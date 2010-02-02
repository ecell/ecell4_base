#ifndef SHELL_HPP
#define SHELL_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "Sphere.hpp"

template<typename T_, typename Tdid_>
struct Shell
{
    typedef Sphere<T_> sphere_type;
    typedef Tdid_ domain_id_type;
    typedef typename sphere_type::position_type position_type;
    typedef typename sphere_type::length_type length_type;

    Shell(): sphere_(), domain_id_() {}

    Shell(domain_id_type const& domain_id, sphere_type const& sphere)
        : sphere_(sphere), domain_id_(domain_id) {}

    length_type calculateDistanceToSelf(position_type pos)
    {
        return sphere_.calculateDistanceToSelf(pos);
    }

    length_type calculateDistanceToSelfWithOffset(position_type pos, 
                                                  position_type offset)
    {
        return sphere_.calculateDistanceToSelfWithOffset(pos, offset);
    }

    position_type& position()
    {
        return sphere_.position();
    }

    position_type const& position() const
    {
        return sphere_.position();
    }

    length_type& radius()
    {
        return sphere_.radius();
    }

    length_type const& radius() const
    {
        return sphere_.radius();
    }

    sphere_type& as_sphere()
    {
        return sphere_;
    }

    sphere_type const& as_sphere() const
    {
        return sphere_;
    }

    domain_id_type const& did() const
    {
        return domain_id_;
    }

    domain_id_type& did()
    {
        return domain_id_;
    }

    bool operator==(Shell const& rhs) const
    {
        return domain_id_ == rhs.did() && sphere_ == rhs.as_sphere();
    }

    bool operator!=(Shell const& rhs) const
    {
        return !operator==(rhs);
    }

private:
    sphere_type sphere_;
    domain_id_type domain_id_;
};

template<typename Tstrm_, typename Ttraits_, typename T_, typename Tdid_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Shell<T_, Tdid_>& v)
{
    strm << "Shell(" << v.as_sphere() << ", " << v.did() << ")";
    return strm;
}

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_, typename Tdid_>
struct hash<Shell<T_, Tdid_> >
{
    typedef Shell<T_, Tdid_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
            hash<typename argument_type::domain_id_type>()(val.did());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* SHELL_HPP */
