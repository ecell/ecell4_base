#ifndef __ECELL4__VOXEL_HPP
#define __ECELL4__VOXEL_HPP

#include <vector>
#include "Identifier.hpp"

namespace ecell4
{

class MolecularType;

struct Voxel
{
    Voxel(ParticleID id, Integer coord, MolecularType* ptr_mt) :
        id(id),
        coord(coord),
        ptr_mt(ptr_mt),
        diffuse_size(0)
    {
    }

    ParticleID id;
    Integer coord;
    Integer diffuse_size;
    std::vector<Voxel*> adjoiningVoxels;
    MolecularType* ptr_mt;

};


} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

