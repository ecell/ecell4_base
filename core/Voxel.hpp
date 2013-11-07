#ifndef __ECELL4__VOXEL_HPP
#define __ECELL4__VOXEL_HPP

#include <vector>
#include "Identifier.hpp"

namespace ecell4
{

class MolecularTypeBase;

struct Voxel
{

public:
    Voxel(ParticleID id, Integer coord, MolecularTypeBase* ptr_mt) :
        id(id),
        coord(coord),
        ptr_mt(ptr_mt)
    {
    }
    ~Voxel()
    {
    }

public:
    ParticleID id;
    Integer coord;
    MolecularTypeBase* ptr_mt;

};

} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

