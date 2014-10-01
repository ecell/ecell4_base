#ifndef __ECELL4_EGFRD_PARTICLE_TRAITS_HPP
#define __ECELL4_EGFRD_PARTICLE_TRAITS_HPP

#include <ecell4/core/Particle.hpp>
#include "Sphere.hpp"
#include "Shape.hpp"

inline
Sphere shape(ecell4::Particle &p)
{
    return Sphere(p.position(), p.radius());
}

inline
Sphere shape(const ecell4::Particle &p)
{
    return Sphere(p.position(), p.radius());
}

#endif
