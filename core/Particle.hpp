#ifndef __ECELL4_PARTICLE_HPP
#define __ECELL4_PARTICLE_HPP

#include <map>

#include "types.hpp"
#include "Position3.hpp"
#include "Species.hpp"
#include "Identifier.hpp"


namespace ecell4
{

class Particle
{
public:

    Particle()
    {
        ;
    }

    Particle(
        const Species& sp, const Position3& pos, const Real& radius,
        const Real& D)
        : species_(sp), position_(pos), radius_(radius), D_(D)
    {
        ;
    }

    Position3& position()
    {
        return position_;
    }

    const Position3& position() const
    {
        return position_;
    }

    Real& radius()
    {
        return radius_;
    }

    const Real& radius() const
    {
        return radius_;
    }

    Real& D()
    {
        return D_;
    }

    const Real& D() const
    {
        return D_;
    }

    Species& species()
    {
        return species_;
    }

    const Species& species() const
    {
        return species_;
    }

private:

    Species species_;
    Position3 position_;
    Real radius_, D_;
};

} // ecell4

#endif /* __ECELL4_PARTICLE_HPP */
