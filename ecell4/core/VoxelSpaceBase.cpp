#include "VoxelSpaceBase.hpp"
#include "Context.hpp"

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

Integer VoxelSpaceBase::num_molecules(const Species& sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            const boost::shared_ptr<VoxelPool>& vp((*itr).second);
            count += count_voxels(vp) * cnt;
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            const boost::shared_ptr<MoleculePool>& vp((*itr).second);
            count += vp->size() * cnt;
        }
    }
    return count;
}


/*
 * ParticleSpace Traits
 */

std::vector<std::pair<ParticleID, Particle> >
VoxelSpaceBase::list_particles() const
{
    const std::vector<std::pair<ParticleID, ParticleVoxel> > voxels(list_voxels());

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(voxels.size());
    for (std::vector<std::pair<ParticleID, ParticleVoxel> >::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
VoxelSpaceBase::list_particles(const Species& sp) const
{
    const std::vector<std::pair<ParticleID, ParticleVoxel> > voxels(list_voxels(sp));

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(voxels.size());
    for (std::vector<std::pair<ParticleID, ParticleVoxel> >::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
VoxelSpaceBase::list_particles_exact(const Species& sp) const
{
    const std::vector<std::pair<ParticleID, ParticleVoxel> >
        voxels(list_voxels_exact(sp));

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(voxels.size());
    for (std::vector<std::pair<ParticleID, ParticleVoxel> >::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}


/*
 * VoxelSpace Traits
 */

bool VoxelSpaceBase::has_voxel(const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        if (vp->find(pid) != vp->end())
            return true;
    }
    return false;
}

Integer VoxelSpaceBase::num_voxels_exact(const Species& sp) const
{
    {
        voxel_pool_map_type::const_iterator itr(voxel_pools_.find(sp));
        if (itr != voxel_pools_.end())
        {
            return count_voxels((*itr).second);
        }
    }

    {
        molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
        if (itr != molecule_pools_.end())
        {
            const boost::shared_ptr<MoleculePool>& vp((*itr).second);
            return vp->size();  // upcast
        }
    }

    return 0;
}

Integer VoxelSpaceBase::num_voxels(const Species& sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        if (sexp.match((*itr).first))
        {
            const boost::shared_ptr<VoxelPool>& vp((*itr).second);
            count += count_voxels(vp);
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if (sexp.match((*itr).first))
        {
            const boost::shared_ptr<MoleculePool>& vp((*itr).second);
            count += vp->size();
        }
    }
    return count;
}

Integer VoxelSpaceBase::num_voxels() const
{
    Integer count(0);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        const boost::shared_ptr<VoxelPool>& vp((*itr).second);
        count += count_voxels(vp);
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        count += vp->size();
    }

    return count;
}

void VoxelSpaceBase::push_voxels(std::vector<std::pair<ParticleID, ParticleVoxel> >& voxels,
        const boost::shared_ptr<MoleculePool>& voxel_pool,
        const Species& species) const
{
    const std::string location_serial(get_location_serial(voxel_pool));
    for (MoleculePool::const_iterator i(voxel_pool->begin()); i != voxel_pool->end(); ++i)
        voxels.push_back(
                std::make_pair(
                    (*i).pid,
                    ParticleVoxel(species, (*i).coordinate, voxel_pool->radius(),
                        voxel_pool->D(), location_serial)));
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
VoxelSpaceBase::list_voxels() const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        push_voxels(retval, vp, vp->species());
    }

    return retval;
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
VoxelSpaceBase::list_voxels(const Species& sp) const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;
    SpeciesExpressionMatcher sexp(sp);

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
            itr != molecule_pools_.end(); ++itr)
        if (sexp.match((*itr).first))
            push_voxels(retval, (*itr).second, sp);

    return retval;
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
VoxelSpaceBase::list_voxels_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;

    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
        push_voxels(retval, (*itr).second, sp);
    return retval;
}

std::pair<ParticleID, ParticleVoxel>
VoxelSpaceBase::get_voxel(const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        MoleculePool::container_type::const_iterator j(vp->find(pid));
        if (j != vp->end())
            return std::make_pair(
                    pid,
                    ParticleVoxel((*itr).first, (*j).coordinate, vp->radius(),
                        vp->D(), get_location_serial(vp)));
    }
    throw NotFound("voxel not found.");
}

boost::shared_ptr<VoxelPool> VoxelSpaceBase::find_voxel_pool(const Species& sp)
{
    voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second;
    }
    return find_molecule_pool(sp);  // upcast
}

boost::shared_ptr<const VoxelPool> VoxelSpaceBase::find_voxel_pool(const Species& sp) const
{
    voxel_pool_map_type::const_iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second;
    }
    return find_molecule_pool(sp);  // upcast
}

bool VoxelSpaceBase::has_molecule_pool(const Species& sp) const
{
    return (molecule_pools_.find(sp) != molecule_pools_.end());
}

boost::shared_ptr<MoleculePool> VoxelSpaceBase::find_molecule_pool(const Species& sp)
{
    molecule_pool_map_type::iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second;  // upcast
    }
    throw NotFound("MoleculePool not found.");
}

boost::shared_ptr<const MoleculePool> VoxelSpaceBase::find_molecule_pool(const Species& sp) const
{
    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second;  // upcast
    }
    throw NotFound("MoleculePool not found.");
}

} // ecell4
