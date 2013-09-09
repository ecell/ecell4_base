#ifndef __ECELL4_MOLECULE_TYPE_HPP
#define __ECELL4_MOLECULE_TYPE_HPP

#include "type.hpp"
#include "Species.hpp"

namespace ecell4
{

struct SParticle
{
    Integer coord;
    Species& species;
}

}

#endif

