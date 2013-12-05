#include "MolecularType.hpp"

namespace ecell4
{

void MolecularType::addVoxel(Integer coord, ParticleID pid)
{
    container_type::iterator itr(find(coord));
    if (itr != voxels_.end())
    {
        voxels_.erase(itr);
    }
    voxels_.push_back(container_type::value_type(
                std::pair<Integer, ParticleID>(coord, pid)));
}

bool MolecularType::removeVoxel(Integer coord)
{
    container_type::iterator itr(find(coord));
    if (itr != voxels_.end())
    {
        voxels_.erase(itr);
        return true;
    }
    return false;
}

const Species& MolecularType::species() const
{
    return this->species_;
}

const MolecularType::container_type& MolecularType::voxels() const
{
    return voxels_;
}

std::vector<SParticle> MolecularType::sparticles() const
{
    std::vector<SParticle> retval;
    for (container_type::const_iterator itr(begin());
            itr != end(); ++itr)
    {
        retval.push_back(SParticle((*itr).first, &species_));
    }
    return retval;
}

} // ecell4
