#ifndef CYLINDRICAL_SURFACE_HPP
#define CYLINDRICAL_SURFACE_HPP

#include <boost/bind.hpp>
#include "Surface.hpp"
#include "Cylinder.hpp"

template<typename Ttraits_>
class CylindricalSurface
    : public BasicSurfaceImpl<Ttraits_, Cylinder>
{
public:
    typedef BasicSurfaceImpl<Ttraits_, Cylinder> base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::shape_type shape_type;
    typedef typename base_type::rng_type rng_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;

    virtual position_type random_position(rng_type& rng) const
    {
        return ::random_position(base_type::shape(), boost::bind(&rng_type::uniform, rng, -1., 1.));
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return multiply(base_type::shape().unit_z(),
                (rng.uniform_int(0, 1) * 2 - 1) * r);
    }

    virtual position_type bd_displacement(length_type const& sigma, rng_type& rng) const
    {
        return multiply(base_type::shape().unit_z(), rng.gaussian(sigma));
    }

    virtual length_type minimal_distance(length_type const& radius) const
    {
        length_type cylinder_radius = base_type::shape().radius();
        // Return minimal distance *to* surface.
        return (cylinder_radius + radius) * traits_type::MINIMAL_SEPARATION_FACTOR - cylinder_radius;
    }

    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    CylindricalSurface(identifier_type const& id, shape_type const& shape)
        : base_type(id, shape) {}
};

#endif /* CYLINDRICAL_SURFACE_HPP */
