#ifndef __ECELL4__VOXEL_HPP
#define __ECELL4__VOXEL_HPP

#include <vector>
#include "Identifier.hpp"

namespace ecell4
{

class MolecularType;

struct Voxel
{
public:
    Voxel(ParticleID id, Integer coord, MolecularType* ptr_mt) :
        id(id),
        coord(coord),
        ptr_mt(ptr_mt),
        diffuse_size(0),
        adjoining_size(12),
        adjoiningVoxels(adjoining_size)
    {
    }
    ~Voxel()
    {
    }
    bool setAdjoiningVoxel(Integer direction, Voxel* adjoining)
    {
        if (0 > direction || direction >= adjoining_size)
        {
            return false;
        }
        adjoiningVoxels[direction] = adjoining;
        return true;
    }

public:
    ParticleID id;
    Integer coord;
    Integer diffuse_size;
    Integer adjoining_size;
    std::vector<Voxel*> adjoiningVoxels;
    MolecularType* ptr_mt;
};


} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

