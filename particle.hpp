#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "species_id.hpp"
#include "sphere.hpp"

template<typename T_>
struct particle
{
    typedef sphere<T_> sphere_type;
    typedef species_id species_id_type;
    typedef typename sphere_type::position_type position_type;
    typedef typename sphere_type::length_type length_type;

    particle(): sphere_(), species_id_() {}

    particle(species_id_type const& species_id, sphere_type const& sphere)
        : sphere_(sphere), species_id_(species_id) {}

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

    species_id_type const& sid() const
    {
        return species_id_;
    }

    species_id_type& sid()
    {
        return species_id_;
    }

private:
    sphere_type sphere_;
    species_id species_id_;
};

#endif /* PARTICLE_HPP */
