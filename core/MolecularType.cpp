#include "MolecularType.hpp"

namespace ecell4
{

void MolecularType::addVoxel(particle_info info)
{
    container_type::iterator itr(find(info.first));
    if (itr != voxels_.end())
    {
        voxels_.erase(itr);
    }
    voxels_.push_back(info);
}

bool MolecularType::removeVoxel(Coord coord)
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
