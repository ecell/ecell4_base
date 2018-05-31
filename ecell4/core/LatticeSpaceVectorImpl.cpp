#include "Context.hpp"
#include "MoleculePool.hpp"
#include "VacantType.hpp"
#include "StructureType.hpp"
#include "LatticeSpaceVectorImpl.hpp"

namespace ecell4 {

typedef LatticeSpaceVectorImpl::coordinate_type coordinate_type;

LatticeSpaceVectorImpl::LatticeSpaceVectorImpl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const bool is_periodic) :
    base_type(edge_lengths, voxel_radius, is_periodic), is_periodic_(is_periodic)
{
    std::stringstream ss;
    ss << voxel_radius_;
    border_ = boost::shared_ptr<VoxelPool>(
            new MoleculePool(Species("Border", ss.str(), "0"), vacant_));
    periodic_ = boost::shared_ptr<VoxelPool>(
            new MoleculePool(Species("Periodic", ss.str(), "0"), vacant_));

    initialize_voxels(is_periodic_);
}

LatticeSpaceVectorImpl::~LatticeSpaceVectorImpl() {}

void LatticeSpaceVectorImpl::initialize_voxels(const bool is_periodic)
{
    const coordinate_type voxel_size(col_size_ * row_size_ * layer_size_);
    // std::cout << "voxel_size = " << voxel_size << std::endl;

    voxel_pools_.clear();
    molecule_pools_.clear();
    voxels_.clear();
    voxels_.reserve(voxel_size);
    for (coordinate_type coord(0); coord < voxel_size; ++coord)
    {
        if (!is_inside(coord))
        {
            if (is_periodic)
            {
                voxels_.push_back(periodic_);
                periodic_->add_voxel(coordinate_id_pair_type(ParticleID(), coord));
            }
            else
            {
                voxels_.push_back(border_);
                border_->add_voxel(coordinate_id_pair_type(ParticleID(), coord));
            }
        }
        else
        {
            voxels_.push_back(vacant_);
            vacant_->add_voxel(coordinate_id_pair_type(ParticleID(), coord));
        }
    }
}

Integer LatticeSpaceVectorImpl::num_species() const
{
    return voxel_pools_.size() + molecule_pools_.size();
}

std::pair<ParticleID, ParticleVoxel>
LatticeSpaceVectorImpl::get_voxel_at(const coordinate_type& coord) const
{
    boost::shared_ptr<const VoxelPool> vp(voxels_.at(coord));

    return std::make_pair(
        vp->get_particle_id(coord),
        ParticleVoxel(vp->species(),
              coord,
              vp->radius(),
              vp->D(),
              get_location_serial(vp)));
}

bool LatticeSpaceVectorImpl::update_structure(const Particle& p)
{
    //XXX: Particle does not have a location.
    ParticleVoxel v(p.species(), position2coordinate(p.position()), p.radius(), p.D());
    return update_voxel(ParticleID(), v);
}

/*
 * original methods
 */

const Species& LatticeSpaceVectorImpl::find_species(std::string name) const
{
    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        if ((*itr).first.serial() == name)
        {
            return (*itr).first;
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if ((*itr).first.serial() == name)
        {
            return (*itr).first;
        }
    }
    throw NotFound(name);
}

std::vector<coordinate_type>
LatticeSpaceVectorImpl::list_coords_exact(const Species& sp) const
{
    std::vector<coordinate_type> retval;

    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr == molecule_pools_.end())
    {
        return retval;
    }

    const boost::shared_ptr<MoleculePool>& vp((*itr).second);

    for (MoleculePool::const_iterator itr(vp->begin()); itr != vp->end(); ++itr)
    {
        retval.push_back((*itr).coordinate);
    }
    return retval;
}

std::vector<coordinate_type>
LatticeSpaceVectorImpl::list_coords(const Species& sp) const
{
    std::vector<coordinate_type> retval;
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if (!spmatch(sp, (*itr).first))
        {
            continue;
        }

        const boost::shared_ptr<MoleculePool>& vp((*itr).second);

        for (MoleculePool::const_iterator vitr(vp->begin());
             vitr != vp->end(); ++vitr)
        {
            retval.push_back((*vitr).coordinate);
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
LatticeSpaceVectorImpl::list_voxels() const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);

        const std::string loc(get_location_serial(vp));
        const Species& sp(vp->species());

        for (MoleculePool::const_iterator i(vp->begin());
            i != vp->end(); ++i)
        {
            retval.push_back(std::make_pair(
                (*i).pid,
                ParticleVoxel(sp, (*i).coordinate, vp->radius(), vp->D(), loc)));
        }
    }

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        const boost::shared_ptr<VoxelPool>& vp((*itr).second);

        const std::string loc(get_location_serial(vp));
        const Species& sp(vp->species());

        for (voxel_container::const_iterator i(voxels_.begin()); i != voxels_.end(); ++i)
        {
            if (*i != vp)
            {
                continue;
            }

            const coordinate_type coord(std::distance(voxels_.begin(), i));
            retval.push_back(std::make_pair(
                ParticleID(),
                ParticleVoxel(sp, coord, vp->radius(), vp->D(), loc)));
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
LatticeSpaceVectorImpl::list_voxels_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;

    {
        voxel_pool_map_type::const_iterator itr(voxel_pools_.find(sp));
        if (itr != voxel_pools_.end())
        {
            const boost::shared_ptr<VoxelPool>& vp((*itr).second);
            const std::string loc(get_location_serial(vp));
            for (voxel_container::const_iterator i(voxels_.begin()); i != voxels_.end(); ++i)
            {
                if (*i != vp)
                {
                    continue;
                }

                const coordinate_type coord(std::distance(voxels_.begin(), i));
                retval.push_back(std::make_pair(
                    ParticleID(),
                    ParticleVoxel(sp, coord, vp->radius(), vp->D(), loc)));
            }
            return retval;
        }
    }

    {
        molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
        if (itr != molecule_pools_.end())
        {
            const boost::shared_ptr<MoleculePool>& vp((*itr).second);
            const std::string loc(get_location_serial(vp));
            for (MoleculePool::const_iterator i(vp->begin());
                 i != vp->end(); ++i)
            {
                retval.push_back(std::make_pair(
                    (*i).pid,
                    ParticleVoxel(sp, (*i).coordinate, vp->radius(), vp->D(), loc)));
            }
            return retval;
        }
    }
    return retval; // an empty vector
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
LatticeSpaceVectorImpl::list_voxels(const Species& sp) const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;
    SpeciesExpressionMatcher sexp(sp);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        if (!sexp.match((*itr).first))
        {
            continue;
        }

        const boost::shared_ptr<VoxelPool>& vp((*itr).second);
        const std::string loc(get_location_serial(vp));
        for (voxel_container::const_iterator i(voxels_.begin()); i != voxels_.end(); ++i)
        {
            if (*i != vp)
            {
                continue;
            }

            const coordinate_type coord(std::distance(voxels_.begin(), i));
            retval.push_back(std::make_pair(
                ParticleID(),
                ParticleVoxel(sp, coord, vp->radius(), vp->D(), loc)));
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if (!sexp.match((*itr).first))
        {
            continue;
        }

        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        const std::string loc(get_location_serial(vp));
        for (MoleculePool::const_iterator i(vp->begin());
            i != vp->end(); ++i)
        {
            retval.push_back(std::make_pair(
                (*i).pid,
                ParticleVoxel(sp, (*i).coordinate, vp->radius(), vp->D(), loc)));
        }
    }

    return retval;
}

/*
 * Protected functions
 */

coordinate_type LatticeSpaceVectorImpl::get_coord(
    const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        for (MoleculePool::const_iterator vitr(vp->begin());
             vitr != vp->end(); ++vitr)
        {
            if ((*vitr).pid == pid)
            {
                return (*vitr).coordinate;
            }
        }
    }
    return -1; //XXX: a bit dirty way
}

// bool LatticeSpaceVectorImpl::has_species_exact(const Species& sp) const
// {
//     return spmap_.find(sp) != spmap_.end();
// }

bool LatticeSpaceVectorImpl::remove_voxel(const ParticleID& pid)
{
    for (molecule_pool_map_type::iterator i(molecule_pools_.begin());
         i != molecule_pools_.end(); ++i)
    {
        const boost::shared_ptr<MoleculePool>& vp((*i).second);
        MoleculePool::const_iterator j(vp->find(pid));
        if (j != vp->end())
        {
            const coordinate_type coord((*j).coordinate);
            if (!vp->remove_voxel_if_exists(coord))
            {
                return false;
            }

            voxels_.at(coord) = vp->location();

            vp->location()->add_voxel(
                coordinate_id_pair_type(ParticleID(), coord));
            return true;
        }
    }
    return false;
}

bool LatticeSpaceVectorImpl::remove_voxel(const coordinate_type& coord)
{
    boost::shared_ptr<VoxelPool> vp(voxels_.at(coord));
    if (vp->is_vacant())
    {
        return false;
    }
    if (vp->remove_voxel_if_exists(coord))
    {
        voxels_.at(coord) = vp->location();
        vp->location()->add_voxel(
            coordinate_id_pair_type(ParticleID(), coord));
        return true;
    }
    return false;
}

bool
LatticeSpaceVectorImpl::move(
        const coordinate_type& src,
        const coordinate_type& dest,
        const std::size_t candidate)
{
    return move_(src, dest, candidate).second;
}

bool LatticeSpaceVectorImpl::can_move(
    const coordinate_type& src, const coordinate_type& dest) const
{
    if (src == dest)
        return false;

    boost::shared_ptr<const VoxelPool> src_vp(voxels_.at(src));
    if (src_vp->is_vacant())
        return false;

    boost::shared_ptr<VoxelPool> dest_vp(voxels_.at(dest));

    if (dest_vp == border_)
        return false;

    if (dest_vp == periodic_)
        dest_vp = voxels_.at(apply_boundary_(dest));

    return (dest_vp == src_vp->location());
}

std::pair<coordinate_type, bool>
LatticeSpaceVectorImpl::move_(
        coordinate_type from,
        coordinate_type to,
        const std::size_t candidate)
{
    if (from == to)
    {
        return std::pair<coordinate_type, bool>(from, false);
    }

    boost::shared_ptr<VoxelPool> from_vp(voxels_.at(from));
    if (from_vp->is_vacant())
    {
        return std::pair<coordinate_type, bool>(from, true);
    }

    boost::shared_ptr<VoxelPool> to_vp(voxels_.at(to));

    if (to_vp == border_)
    {
        return std::pair<coordinate_type, bool>(from, false);
    }
    else if (to_vp == periodic_)
    {
        to = apply_boundary_(to);
        to_vp = voxels_.at(to);
    }

    if (to_vp != from_vp->location())
    {
        return std::pair<coordinate_type, bool>(to, false);
    }

    from_vp->replace_voxel(from, to, candidate);
    voxels_.at(from) = to_vp;

    to_vp->replace_voxel(to, from);
    voxels_.at(to) = from_vp;

    return std::pair<coordinate_type, bool>(to, true);
}

std::pair<coordinate_type, bool>
LatticeSpaceVectorImpl::move_(
        coordinate_id_pair_type& info,
        coordinate_type to)
{
    const coordinate_type from(info.coordinate);
    if (from == to)
    {
        return std::pair<coordinate_type, bool>(from, false);
    }

    boost::shared_ptr<VoxelPool> from_vp(voxels_.at(from));
    if (from_vp->is_vacant())
    {
        return std::pair<coordinate_type, bool>(from, true);
    }

    boost::shared_ptr<VoxelPool> to_vp(voxels_.at(to));

    if (to_vp == border_)
    {
        return std::pair<coordinate_type, bool>(from, false);
    }
    else if (to_vp == periodic_)
    {
        to = apply_boundary_(to);
        to_vp = voxels_.at(to);
    }

    if (to_vp != from_vp->location())
    {
        return std::pair<coordinate_type, bool>(to, false);
    }

    info.coordinate = to;
    voxels_.at(from) = to_vp;

    // to_vp->replace_voxel(to, coordinate_id_pair_type(ParticleID(), from));
    to_vp->replace_voxel(to, from);
    voxels_.at(to) = from_vp;

    return std::pair<coordinate_type, bool>(to, true);
}

const Particle LatticeSpaceVectorImpl::particle_at(
    const coordinate_type& coord) const
{
    boost::shared_ptr<const VoxelPool> vp(voxels_.at(coord));

    return Particle(vp->species(),
                    coordinate2position(coord),
                    vp->radius(),
                    vp->D());
}

/*
 * Change the Species and coordinate of a ParticleVoxel with ParticleID, pid, to
 * v.species() and v.coordinate() respectively and return false.
 * If no ParticleVoxel with pid is found, create a new ParticleVoxel at v.coordiante() and return ture.
 */
bool LatticeSpaceVectorImpl::update_voxel(const ParticleID& pid, ParticleVoxel v)
{
    const coordinate_type& to_coord(v.coordinate);
    if (!is_in_range(to_coord))
    {
        throw NotSupported("Out of bounds");
    }

    boost::shared_ptr<VoxelPool> new_vp(get_voxel_pool(v)); //XXX: need MoleculeInfo
    boost::shared_ptr<VoxelPool> dest_vp(get_voxel_pool_at(to_coord));

    if (dest_vp != new_vp->location())
    {
        throw NotSupported(
            "Mismatch in the location. Failed to place '"
            + new_vp->species().serial() + "' to '"
            + dest_vp->species().serial() + "'.");
    }

    const coordinate_type
        from_coord(pid != ParticleID() ? get_coord(pid) : -1);
    if (from_coord != -1)
    {
        // move
        voxels_.at(from_coord)
               ->remove_voxel_if_exists(from_coord);

        //XXX: use location?
        dest_vp->replace_voxel(to_coord, from_coord);
        voxels_.at(from_coord) = dest_vp;

        new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));
        voxels_.at(to_coord) = new_vp;
        return false;
    }

    // new
    dest_vp->remove_voxel_if_exists(to_coord);

    new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));
    voxels_.at(to_coord) = new_vp;
    return true;
}

bool
LatticeSpaceVectorImpl::add_voxel(
        const Species& sp,
        const ParticleID& pid,
        const coordinate_type& coordinate)
{
    boost::shared_ptr<VoxelPool> vpool(find_voxel_pool(sp));
    boost::shared_ptr<VoxelPool> location(get_voxel_pool_at(coordinate));

    if (vpool->location() != location)
        return false;

    location->remove_voxel_if_exists(coordinate);
    vpool->add_voxel(coordinate_id_pair_type(pid, coordinate));
    voxels_.at(coordinate) = vpool;

    return true;
}

bool
LatticeSpaceVectorImpl::add_voxels(
        const Species& sp,
        std::vector<std::pair<ParticleID, coordinate_type> > voxels)
{
    // this function doesn't check location.
    boost::shared_ptr<VoxelPool> mtb;
    try
    {
        mtb = find_voxel_pool(sp);
    }
    catch (NotFound &e)
    {
        return false;
    }

    for (std::vector<std::pair<ParticleID, coordinate_type> >::iterator itr(voxels.begin());
            itr != voxels.end(); ++itr)
    {
        const ParticleID pid((*itr).first);
        const coordinate_type coord((*itr).second);
        get_voxel_pool_at(coord)->remove_voxel_if_exists(coord);
        mtb->add_voxel(coordinate_id_pair_type(pid, coord));
        voxels_.at(coord) = mtb;
    }
    return true;
}

} // ecell4
