#ifndef CYLINDRICAL_SURFACE_HPP
#define CYLINDRICAL_SURFACE_HPP

#include "Surface.hpp"
#include "Cylinder.hpp"

template<typename Ttraits_>
class CylindricalSurface
    : public BasicRegionImpl<Ttraits_, Cylinder<typename Ttraits_::length_type> >
{
public:
    typedef BasicRegionImpl<Ttraits_, Cylinder<typename Ttraits_::length_type> > base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::shape_type shape_type;
    typedef typename base_type::rng_type rng_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;

    virtual position_type random_position(rng_type& rng) const
    {
        return add(
            base_type::shape().position(),
            base_type::shape().unit_z() * rng.uniform(-1., 1.));
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return multiply(base_type::shape().unit_z(),
                (rng.uniform_int(0, 1) * 2 - 1) * r);
    }

    virtual position_type bd_displacement(length_type const& r, rng_type& rng) const
    {
        return multiply(base_type::shape().unit_z(), rng.normal(0., r));
    }

    virtual length_type minimal_distance(length_type const& radius) const
    {
        return (base_type::shape().radius() + radius) * traits_type::MINIMAL_SEPARATION_FACTOR;
    }

    CylindricalSurface(identifier_type const& id, shape_type const& shape)
        : base_type(id, shape) {}
};

#endif /* CYLINDRICAL_SURFACE_HPP */
