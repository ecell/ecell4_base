#include "MolecularType.hpp"

namespace ecell4
{

void MolecularType::addVoxel(Voxel &voxel)
{
    voxels_.insert(voxel_container_type::value_type(
                voxel.id, voxel));
    voxel.ptr_mt = this;
}

bool MolecularType::removeVoxel(const ParticleID pid)
{
    voxel_container_type::iterator itr(voxels_.find(pid));
    if (itr == voxels_.end())
    {
        return false;
    }
    (*itr).second.ptr_mt = NULL;
    voxels_.erase(itr++);
    return true;
}

const Species& MolecularType::species() const
{
    return this->species_;
}

const MolecularType::voxel_container_type& MolecularType::voxels() const
{
    return voxels_;
}

}
