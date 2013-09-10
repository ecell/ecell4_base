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
    std::vector<Voxel> adjoiningVoxels;
    MolecularType *p_molecule_type;

};
bool compare_Voxel(const Voxel &v1, const Voxel &v2)
{
    if (v1.id == v2.id) {
        if (v1.coord == v2.coord) {
            if (v1.diffuseSize == v2.diffuseSize) {
                return true;
            }
        }
    }
    return false;
}

} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

