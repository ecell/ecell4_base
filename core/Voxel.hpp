#ifndef __ECELL4__VOXEL_HPP
#define __ECELL4__VOXEL_HPP

#include <vector>
#include "Identifier.hpp"

namespace ecell4
{

class MolecularType;

struct Voxel
{

    ParticleID id;
    Integer coord;
    Integer diffuse_size;
    std::vector<Voxel*> adjoiningVoxels;
    MolecularType *p_molecule_type;

};


} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

