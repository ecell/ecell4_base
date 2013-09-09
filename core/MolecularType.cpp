#include "MolecularType.hpp"

namespace ecell4
{

MolecularType::MolecularType(Species& species)
{
    species_ = species;
}

void MolecularType::addVoxel(const Voxel& voxel)
{
    voxels_.push_back(voxel);
    voxel.molecule_type = static_cast<MolecularType&>(*this);
}

void MolecularType::removeVoxel(const Voxel& voxel)
{
    for (voxel_container_type::iterator i(voxels_.begin());
            i != voxels_.end(); ++i)
    {
        if (*i == voxel)
            voxels_.erase(i);
    }
}

const Species& MolecularType::species() const
{
    return sp_;
}

const voxel_container_type MolecularType::voxels() const
{
    return voxels_;
}

}
