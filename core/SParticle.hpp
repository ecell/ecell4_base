#ifndef __ECELL4_SPARTICLE_HPP
#define __ECELL4_SPARTICLE_HPP

#include "types.hpp"
#include "Species.hpp"

namespace ecell4
{

struct SParticle
{
    Integer coord;
    Species& species;
};

}

#endif

