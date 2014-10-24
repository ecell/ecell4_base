#ifndef __ECELL4_EGFRD_PARTICLE_TRAITS_HPP
#define __ECELL4_EGFRD_PARTICLE_TRAITS_HPP

#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Sphere.hpp>
#include "Sphere.hpp"
#include "Shape.hpp"

inline
ecell4::Sphere shape(ecell4::Particle &p)
{
    return ecell4::Sphere(p.position(), p.radius());
    // return Sphere(p.position(), p.radius());
}

inline
ecell4::Sphere shape(const ecell4::Particle &p)
{
    return ecell4::Sphere(p.position(), p.radius());
    // return Sphere(p.position(), p.radius());
}

inline ecell4::Sphere::length_type
distance(const ecell4::Sphere& obj, const ecell4::Sphere::position_type& pos)
{
    return distance(pos, obj.position()) - obj.radius();
}

template<typename T_>
inline ecell4::Sphere::length_type
distance_cyclic(
    const ecell4::Sphere& p1, T_ const& p2,
    const ecell4::Sphere::position_type& edge_lengths)
{
    return distance(p1, cyclic_transpose(p2, p1.position(), edge_lengths));
}

template<>
struct shape_position_type<ecell4::Sphere> {
    typedef ecell4::Sphere::position_type type;
};

namespace ecell4
{

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const ecell4::Sphere& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << "}";
    return strm;
}

} // ecell4

#endif
