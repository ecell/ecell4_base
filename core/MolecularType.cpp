#include "MolecularType.hpp"

namespace ecell4
{

void MolecularType::addVoxel(Voxel *voxel, ParticleID pid)
{
    if (find(pid) == voxels_.end())
    {
        voxels_.push_back(container_type::value_type(
                    std::pair<Voxel*, ParticleID>(voxel, pid)));
        voxel->ptr_mt = this;
    }
}

bool MolecularType::removeVoxel(const ParticleID pid)
{
    container_type::iterator itr(find(pid));
    bool flg(false);
    if (itr != voxels_.end())
    {
        voxels_.erase(itr);
        (*itr).first->ptr_mt = NULL;
        flg = true;
    }
    return flg;
}

const Species& MolecularType::species() const
{
    return this->species_;
}

const MolecularType::container_type& MolecularType::voxels() const
{
    return voxels_;
}

} // ecell4
