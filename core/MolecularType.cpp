#include "MolecularType.hpp"

namespace ecell4
{

MolecularType::MolecularType(Species& species)
{
    species_ = species;
}

void MolecularType::addVoxel(Voxel& voxel)
{
    voxels_.push_back(voxel);
    voxel.p_molecule_type = this; //static_cast<MolecularType&>(*this);
}

void MolecularType::removeVoxel(const Voxel& voxel)
{
    for (voxel_container_type::iterator i(this->voxels_.begin());
            i != voxels_.end(); ++i)
    {
        if ( i->id == voxel.id )
            voxels_.erase(i);
    }
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
