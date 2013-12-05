#ifndef __ECELL4_SPARTICLE_HPP
#define __ECELL4_SPARTICLE_HPP

#include "types.hpp"
#include "Species.hpp"

namespace ecell4
{

struct SParticle
{
    Integer coord;
    const Species* species;

    SParticle(Integer coord, const Species* species)
        : coord(coord), species(species)
    {
    }
};

} // ecell4

#endif

