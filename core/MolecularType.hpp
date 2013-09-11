#ifndef __ECELL4_MOLECULE_TYPE_HPP
#define __ECELL4_MOLECULE_TYPE_HPP

#include <vector>
#include "Species.hpp"
#include "Voxel.hpp"

namespace ecell4
{

class MolecularType
{

public:

    typedef std::vector<Voxel*> voxel_container_type;

public:

    MolecularType(const std::string& name = "")
        : species_(name)
    {
    }
    void addVoxel(Voxel* voxel);
    void removeVoxel(const Voxel& voxel);
    const Species& species() const;
    const voxel_container_type& voxels() const;

protected:

    Species species_;
    voxel_container_type voxels_;

};

}

#endif
