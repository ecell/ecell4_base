#ifndef SPHERICAL_SURFACE_HPP
#define SPHERICAL_SURFACE_HPP

#include "Surface.hpp"
#include "Sphere.hpp"

template<typename Ttraits_>
class SphericalSurface
    : public BasicRegionImpl<Ttraits_, Sphere<typename Ttraits_::length_type> >
{
public:
    typedef BasicRegionImpl<Ttraits_, Sphere<typename Ttraits_::length_type> > base_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::shape_type shape_type;
    typedef typename base_type::rng_type rng_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;

    virtual position_type random_position(rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual position_type bd_displacement(length_type const& r, rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual length_type minimal_distance(length_type const& radius) const
    {
        return 0.; // TODO
    }

    SphericalSurface(identifier_type const& id, shape_type const& shape)
        : base_type(id, shape) {}
};

#endif /* SPHERICAL_SURFACE_HPP */
