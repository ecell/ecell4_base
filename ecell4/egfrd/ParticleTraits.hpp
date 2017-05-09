#ifndef ECELL4_EGFRD_PARTICLE_TRAITS_HPP
#define ECELL4_EGFRD_PARTICLE_TRAITS_HPP

#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Sphere.hpp>
#include <ecell4/core/Cylinder.hpp>
// #include "Sphere.hpp"
// #include "Shape.hpp"
//XXX: Shape.hpp
#include <boost/type_traits/remove_cv.hpp>
#include "Real3Type.hpp"
#include "geometry.hpp"
//XXX

inline ecell4::Particle offset(
    ecell4::Particle const& shape, ecell4::Particle::position_type off)
{
    ecell4::Particle retval(shape);
    retval.position() += off;
    return retval;
}

// inline
// ecell4::Sphere shape(ecell4::Particle &p)
// {
//     return ecell4::Sphere(p.position(), p.radius());
// }

inline
ecell4::Sphere shape(const ecell4::Particle &p)
{
    return ecell4::Sphere(p.position(), p.radius());
}

// inline ecell4::Sphere offset(
//     const ecell4::Sphere& shape, ecell4::Sphere::position_type off)
// {
//     ecell4::Sphere retval(shape);
//     retval.position() += off;
//     return retval;
// }

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
    return distance(p1, periodic_transpose(p2, p1.position(), edge_lengths));
}

inline ecell4::Sphere::length_type const& shape_size(ecell4::Sphere const& shape)
{
    return shape.size();
}

inline ecell4::Sphere::length_type& shape_size(ecell4::Sphere &shape)
{
    return shape.size();
}

inline ecell4::Sphere::position_type const& shape_position(ecell4::Sphere const& shape)
{
    return shape.position();
}

inline ecell4::Sphere::position_type& shape_position(ecell4::Sphere &shape)
{
    return shape.position();
}

template<>
struct shape_position_type<ecell4::Sphere> {
    typedef ecell4::Sphere::position_type type;
};

// inline ecell4::Cylinder offset(
//     const ecell4::Cylinder& shape, ecell4::Cylinder::position_type off)
// {
//     ecell4::Cylinder retval(shape);
//     retval.position() += off;
//     return retval;
// }

inline ecell4::Cylinder::length_type
distance(const ecell4::Cylinder& obj, const ecell4::Cylinder::position_type& pos)
{
    return distance(pos, obj.position()) - obj.radius();
}

template<typename T_>
inline ecell4::Cylinder::length_type
distance_cyclic(
    const ecell4::Cylinder& p1, T_ const& p2,
    const ecell4::Cylinder::position_type& edge_lengths)
{
    return distance(p1, periodic_transpose(p2, p1.position(), edge_lengths));
}

inline ecell4::Cylinder::length_type const& shape_size(ecell4::Cylinder const& shape)
{
    return shape.size();
}

inline ecell4::Cylinder::length_type& shape_size(ecell4::Cylinder &shape)
{
    return shape.size();
}

inline ecell4::Cylinder::position_type const& shape_position(ecell4::Cylinder const& shape)
{
    return shape.position();
}

inline ecell4::Cylinder::position_type& shape_position(ecell4::Cylinder &shape)
{
    return shape.position();
}

template<>
struct shape_position_type<ecell4::Cylinder> {
    typedef ecell4::Cylinder::position_type type;
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

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const ecell4::Cylinder& v)
{
    strm << "{" << v.position() <<  ", " << v.radius()
        << ", " << v.axis() << ", " << v.half_height() << "}";
    return strm;
}

} // ecell4

#endif
