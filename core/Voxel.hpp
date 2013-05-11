#ifndef __ECELL4__VOXEL_HPP
#define __ECELL4__VOXEL_HPP

#include "Species.hpp";
#include "Identifier.hpp";


namespace ecell4
{

struct Voxel
{
    ParticleID id;
    Species species;
};

} // ecell4

#endif /* __ECELL4__VOXEL_HPP */
