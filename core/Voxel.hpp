#ifndef __ECELL4__VOXEL_HPP
#define __ECELL4__VOXEL_HPP

#include <vector>
#include "MolecularType.hpp"
#include "Identifier.hpp"

namespace ecell4
{

class MolecularType;

struct Voxel
{

    ParticleID id;
    Integer coord;
    Integer diffuseSize;
    std::vector<Voxel&> adjoiningVoxels;
    MolecularType& molecule_type;

};

} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

