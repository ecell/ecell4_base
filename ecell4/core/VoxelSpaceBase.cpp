#include "Context.hpp"
#include "VoxelSpaceBase.hpp"
#include "get_mapper_mf.hpp" // for retrieve_keys(map, keys)

namespace ecell4
{

/*
 * CompartmentSpace Traits
 */

std::vector<Species> VoxelSpaceBase::list_species() const
{
    std::vector<Species> keys;
    utils::retrieve_keys(voxel_pools_, keys);
    utils::retrieve_keys(molecule_pools_, keys);
    return keys;
}

Integer VoxelSpaceBase::num_molecules(const Species &sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    {
        const auto cnt = sexp.count(vacant_->species());
        if (cnt > 0)
        {
            count += vacant_->size() * cnt;
        }
    }

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            count += itr->second->size() * cnt;
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            const boost::shared_ptr<MoleculePool> &vp((*itr).second);
            count += vp->size() * cnt;
        }
    }
    return count;
}

/*
 * VoxelSpace Traits
 */

bool VoxelSpaceBase::has_voxel(const ParticleID &pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool> &vp((*itr).second);
        if (vp->find(pid) != vp->end())
            return true;
    }
    return false;
}

Integer VoxelSpaceBase::num_voxels_exact(const Species &sp) const
{
    if (sp == vacant_->species())
    {
        return vacant_->size();
    }

    {
        voxel_pool_map_type::const_iterator itr(voxel_pools_.find(sp));
        if (itr != voxel_pools_.end())
        {
            return itr->second->size();
        }
    }

    {
        molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
        if (itr != molecule_pools_.end())
        {
            const boost::shared_ptr<MoleculePool> &vp((*itr).second);
            return vp->size();
        }
    }

    return 0;
}

Integer VoxelSpaceBase::num_voxels(const Species &sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    if (sexp.match(vacant_->species()))
    {
        count += vacant_->size();
    }

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        if (sexp.match((*itr).first))
        {
            count += itr->second->size();
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if (sexp.match((*itr).first))
        {
            const boost::shared_ptr<MoleculePool> &vp((*itr).second);
            count += vp->size();
        }
    }
    return count;
}

Integer VoxelSpaceBase::num_voxels() const
{
    Integer count(0);

    if (vacant_->species().serial() != "")
    {
        count += vacant_->size();
    }

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        count += itr->second->size();
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool> &vp((*itr).second);
        count += vp->size();
    }

    return count;
}

static inline void
push_voxels(std::vector<VoxelView> &voxels,
            const boost::shared_ptr<const MoleculePool> &voxel_pool)
{
    for (const auto &voxel : *voxel_pool)
    {
        voxels.push_back(
            VoxelView(voxel.pid, voxel_pool->species(), voxel.coordinate));
    }
}

std::vector<VoxelView> VoxelSpaceBase::list_voxels() const
{
    std::vector<VoxelView> retval;

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool> &vp((*itr).second);
        push_voxels(retval, vp);
    }

    return retval;
}

std::vector<VoxelView> VoxelSpaceBase::list_voxels(const Species &sp) const
{
    std::vector<VoxelView> retval;
    SpeciesExpressionMatcher sexp(sp);

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
        if (sexp.match((*itr).first))
            push_voxels(retval, (*itr).second);

    return retval;
}

std::vector<VoxelView>
VoxelSpaceBase::list_voxels_exact(const Species &sp) const
{
    std::vector<VoxelView> retval;

    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
        push_voxels(retval, (*itr).second);
    return retval;
}

boost::shared_ptr<VoxelPool> VoxelSpaceBase::find_voxel_pool(const Species &sp)
{
    if (sp == vacant_->species())
        return vacant_;

    voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second;
    }
    return find_molecule_pool(sp); // upcast
}

boost::shared_ptr<const VoxelPool>
VoxelSpaceBase::find_voxel_pool(const Species &sp) const
{
    if (sp == vacant_->species())
        return vacant_;

    voxel_pool_map_type::const_iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second;
    }
    return find_molecule_pool(sp); // upcast
}

bool VoxelSpaceBase::has_molecule_pool(const Species &sp) const
{
    return (molecule_pools_.find(sp) != molecule_pools_.end());
}

boost::shared_ptr<MoleculePool>
VoxelSpaceBase::find_molecule_pool(const Species &sp)
{
    molecule_pool_map_type::iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second; // upcast
    }
    throw NotFound("MoleculePool not found.");
}

boost::shared_ptr<const MoleculePool>
VoxelSpaceBase::find_molecule_pool(const Species &sp) const
{
    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second; // upcast
    }
    throw NotFound("MoleculePool not found.");
}

bool VoxelSpaceBase::make_molecular_type(const Species &sp,
                                         const std::string loc)
{
    molecule_pool_map_type::iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return false;
    }
    else if (voxel_pools_.find(sp) != voxel_pools_.end())
    {
        throw IllegalState("The given species is already assigned to the "
                           "VoxelPool with no voxels.");
    }

    boost::shared_ptr<MoleculePool> vp(
        new MoleculePool(sp, find_voxel_pool(Species(loc))));

    std::pair<molecule_pool_map_type::iterator, bool> retval(
        molecule_pools_.insert(molecule_pool_map_type::value_type(sp, vp)));

    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval.second;
}

bool VoxelSpaceBase::make_structure_type(const Species &sp,
                                         const std::string loc)
{
    voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return false;
    }
    else if (molecule_pools_.find(sp) != molecule_pools_.end())
    {
        throw IllegalState(
            "The given species is already assigned to the MoleculePool.");
    }

    boost::shared_ptr<VoxelPool> vp(
        new StructureType(sp, find_voxel_pool(Species(loc))));
    std::pair<voxel_pool_map_type::iterator, bool> retval(
        voxel_pools_.insert(voxel_pool_map_type::value_type(sp, vp)));
    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval.second;
}

} // namespace ecell4
