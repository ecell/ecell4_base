#ifndef __ECELL4__VOXEL_HPP
#define __ECELL4__VOXEL_HPP

#include <vector:c>;
#include "Identifier.hpp";

namespace ecell4
{

class MolecularType;

struct Voxel
{

    ParticleID id(0);
    Integer coord(0);
    Integer diffuseSize(0);
    std::vector<Voxel&> adjoiningVoxels;
    MolecularType& molecule_type;

};

} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

