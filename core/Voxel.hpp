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
    Voxel(Integer coord, MolecularTypeBase* ptr_mt) :
        coord(coord),
        ptr_mt(ptr_mt)
    {
    }
    ~Voxel()
    {
    }

public:
    Integer coord;
    MolecularTypeBase* ptr_mt;

};

} // ecell4

#endif /* __ECELL4__VOXEL_HPP */

