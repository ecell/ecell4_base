#include "LatticeSpaceCellListImpl.hpp"
#include "StructureType.hpp"


namespace ecell4
{

Integer LatticeSpaceCellListImpl::num_molecules(const Species& sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    // for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
    //      itr != voxel_pools_.end(); ++itr)
    // {
    //     const Integer cnt(sexp.count((*itr).first));
    //     if (cnt > 0)
    //     {
    //         const boost::shared_ptr<VoxelPool>& mt((*itr).second);
    //         count += count_voxels(mt) * cnt;
    //     }
    // }

    for (molecular_type_map_type::const_iterator itr(molecular_types_.begin());
         itr != molecular_types_.end(); ++itr)
    {
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            const boost::shared_ptr<MolecularType>& mt((*itr).second);
            count += mt->size() * cnt;
        }
    }
    return count;
}

/*
 * Change the Species and coordinate of a Voxel with ParticleID, pid, to
 * v.species() and v.coordinate() respectively and return false.
 * If no Voxel with pid is found, create a new Voxel at v.coordiante() and return ture.
 */
bool LatticeSpaceCellListImpl::update_voxel(const ParticleID& pid, const Voxel& v)
{
    const coordinate_type& to_coord(v.coordinate());
    if (!is_in_range(to_coord))
    {
        throw NotSupported("Out of bounds");
    }

    VoxelPool* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    VoxelPool* dest_mt(find_molecular_type(to_coord));

    if (dest_mt != new_mt->location())
    {
        throw NotSupported("Mismatch in the location.");
    }

    if (pid != ParticleID())
    {
        const std::pair<VoxelPool*, coordinate_type>
            target(__get_coordinate(pid));
        const coordinate_type& from_coord(target.second);
        if (from_coord != -1)
        {
            // move
            VoxelPool* src_mt(target.first);
            src_mt->remove_voxel_if_exists(from_coord);

            //XXX: use location?
            dest_mt->replace_voxel(to_coord, from_coord);

            new_mt->add_voxel_without_checking(coordinate_id_pair_type(pid, to_coord));

            if (!dest_mt->is_vacant())
            {
                update_matrix(from_coord, dest_mt);
                update_matrix(to_coord, new_mt);
            }
            else
            {
                update_matrix(from_coord, to_coord, new_mt);
            }
            return true;
        }
    }

    // new
    dest_mt->remove_voxel_if_exists(to_coord);

    new_mt->add_voxel_without_checking(coordinate_id_pair_type(pid, to_coord));
    update_matrix(to_coord, new_mt);
    return true;
}

std::pair<VoxelPool*, LatticeSpaceCellListImpl::coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid)
{
    for (molecular_type_map_type::iterator itr(molecular_types_.begin());
         itr != molecular_types_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        MolecularType::container_type::const_iterator j(mt->find(pid));
        if (j != mt->end())
        {
            return std::pair<VoxelPool*, coordinate_type>(
                mt.get(), (*j).coordinate);
        }
    }
    return std::make_pair<VoxelPool*, coordinate_type>(NULL, -1); //XXX: a bit dirty way
}

std::pair<const VoxelPool*, LatticeSpaceCellListImpl::coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid) const
{
    for (molecular_type_map_type::const_iterator itr(molecular_types_.begin());
         itr != molecular_types_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        MolecularType::container_type::const_iterator j(mt->find(pid));
        if (j != mt->end())
        {
            return std::pair<const VoxelPool*, coordinate_type>(
                mt.get(), (*j).coordinate);
        }
    }
    return std::make_pair<const VoxelPool*, coordinate_type>(NULL, -1); //XXX: a bit dirty way
}

VoxelPool* LatticeSpaceCellListImpl::find_molecular_type(
    const LatticeSpaceCellListImpl::coordinate_type& coord) const
{
    /**
     XXX: This may not work
     */
    if (!is_in_range(coord))
    {
        throw NotSupported("Out of bounds");
    }

    if (!is_inside(coord))
    {
        if (is_periodic_)
        {
            return periodic_;
        }
        else
        {
            return border_;
        }
    }

    // for (spmap::iterator itr(spmap_.begin());
    //     itr != spmap_.end(); ++itr)
    // {
    //     VoxelPool& mt((*itr).second);
    //     if (mt.is_vacant())
    //     {
    //         continue;
    //     }

    //     VoxelPool::container_type::const_iterator j(mt.find(coord));
    //     if (j != mt.end())
    //     {
    //         return (&mt);
    //     }
    // }

    const cell_type& cell(matrix_[coordinate2index(coord)]);
    if (cell.size() == 0)
    {
        return vacant_;
    }

    cell_type::const_iterator i(find_from_cell(coord, cell));
    if (i != cell.end())
    {
        return (*i).first;
    }

    return vacant_;
}

VoxelPool* LatticeSpaceCellListImpl::get_molecular_type(const Voxel& v)
{
    const Species& sp(v.species());

    {
        voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
        if (itr != voxel_pools_.end())
        {
            return (*itr).second.get();
        }
    }

    {
        molecular_type_map_type::iterator itr(molecular_types_.find(sp));
        if (itr != molecular_types_.end())
        {
            return (*itr).second.get();  // upcast
        }
    }

    const bool suc = make_molecular_type(sp, v.radius(), v.D(), v.loc());
    if (!suc)
    {
        throw IllegalState("never reach here");
    }

    molecular_type_map_type::iterator i = molecular_types_.find(sp);
    if (i == molecular_types_.end())
    {
        throw IllegalState("never reach here");
    }
    return (*i).second.get();  // upcast
}

bool LatticeSpaceCellListImpl::make_molecular_type(const Species& sp, Real radius, Real D, const std::string loc)
{
    molecular_type_map_type::iterator itr(molecular_types_.find(sp));
    if (itr != molecular_types_.end())
    {
        return false;
    }
    else if (voxel_pools_.find(sp) != voxel_pools_.end())
    {
        throw IllegalState(
            "The given species is already assigned to the VoxelPool with no voxels.");
    }

    VoxelPool* location;
    if (loc == "")
    {
        location = vacant_;
    }
    else
    {
        const Species locsp(loc);
        try
        {
            location = find_molecular_type(locsp);
        }
        catch (const NotFound& err)
        {
            // XXX: A VoxelPool for the structure (location) must be allocated
            // XXX: before the allocation of a Species on the structure.
            // XXX: The VoxelPool cannot be automatically allocated at the time
            // XXX: because its MoleculeInfo is unknown.
            // XXX: LatticeSpaceVectorImpl::load will raise a problem about this issue.
            // XXX: In this implementation, the VoxelPool for a structure is
            // XXX: created with default arguments.
            boost::shared_ptr<MolecularType>
                locmt(new MolecularType(locsp, vacant_, voxel_radius_, 0));
            std::pair<molecular_type_map_type::iterator, bool>
                locval(molecular_types_.insert(
                    molecular_type_map_type::value_type(locsp, locmt)));
            if (!locval.second)
            {
                throw AlreadyExists(
                    "never reach here. find_molecular_type seems wrong.");
            }
            location = (*locval.first).second.get();
        }
    }

    boost::shared_ptr<MolecularType>
        mt(new MolecularType(sp, location, radius, D));
    std::pair<molecular_type_map_type::iterator, bool>
        retval(molecular_types_.insert(
            molecular_type_map_type::value_type(sp, mt)));
    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval.second;
}

std::pair<LatticeSpaceCellListImpl::coordinate_type, bool>
    LatticeSpaceCellListImpl::move_to_neighbor(
        VoxelPool* const& from_mt, VoxelPool* const& loc,
        LatticeSpaceCellListImpl::coordinate_id_pair_type& info, const Integer nrand)
{
    const coordinate_type from(info.coordinate);
    coordinate_type to(get_neighbor(from, nrand));

    VoxelPool* to_mt(find_molecular_type(to));

    if (to_mt != loc)
    {
        if (to_mt == border_)
        {
            return std::make_pair(from, false);
        }
        else if (to_mt != periodic_)
        {
            return std::make_pair(to, false);
        }

        // to_mt == periodic_
        to = periodic_transpose(to);
        to_mt = find_molecular_type(to);

        if (to_mt != loc)
        {
            return std::make_pair(to, false);
        }
    }

    info.coordinate = to; //XXX: updating data

    to_mt->replace_voxel(to, from);

    if (to_mt != vacant_) // (!to_mt->is_vacant())
    {
        update_matrix(from, to_mt);
        update_matrix(to, from_mt);
    }
    else
    {
        update_matrix(from, to, from_mt);
    }
    return std::make_pair(to, true);
}

bool LatticeSpaceCellListImpl::make_structure_type(
    const Species& sp, Shape::dimension_kind dimension, const std::string loc)
{
    voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return false;
    }
    else if (molecular_types_.find(sp) != molecular_types_.end())
    {
        throw IllegalState(
            "The given species is already assigned to the MolecularType.");
    }

    VoxelPool* location;
    if (loc == "")
    {
        location = vacant_;
    }
    else
    {
        const Species locsp(loc);
        try
        {
            location = find_molecular_type(locsp);
        }
        catch (const NotFound& err)
        {
            // XXX: A VoxelPool for the structure (location) must be allocated
            // XXX: before the allocation of a Species on the structure.
            // XXX: The VoxelPool cannot be automatically allocated at the time
            // XXX: because its MoleculeInfo is unknown.
            // XXX: LatticeSpaceVectorImpl::load will raise a problem about this issue.
            // XXX: In this implementation, the VoxelPool for a structure is
            // XXX: created with default arguments.
            boost::shared_ptr<MolecularType>
                locmt(new MolecularType(locsp, vacant_, voxel_radius_, 0));
            std::pair<molecular_type_map_type::iterator, bool>
                locval(molecular_types_.insert(
                    molecular_type_map_type::value_type(locsp, locmt)));
            if (!locval.second)
            {
                throw AlreadyExists(
                    "never reach here. make_structure_type seems wrong.");
            }
            location = (*locval.first).second.get();
        }
    }

    boost::shared_ptr<VoxelPool>
        mt(new StructureType(sp, location, voxel_radius_, dimension));
    std::pair<voxel_pool_map_type::iterator, bool>
        retval(voxel_pools_.insert(voxel_pool_map_type::value_type(sp, mt)));
    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval.second;
}

void LatticeSpaceCellListImpl::add_structure(const Species& sp,
    const boost::shared_ptr<const Shape>& s, const std::string loc)
{
    make_structure_type(sp, s->dimension(), loc);

    structure_container_type::const_iterator i(structures_.find(sp));
    if (i != structures_.end())
    {
        throw NotSupported("not supported yet.");
    }
    structures_.insert(std::make_pair(sp, s));
}

const boost::shared_ptr<const Shape>& LatticeSpaceCellListImpl::get_structure(
    const Species& sp) const
{
    structure_container_type::const_iterator i(structures_.find(sp));
    if (i == structures_.end())
    {
        throw NotFound("not found.");
    }
    return (*i).second;
}

const Shape::dimension_kind LatticeSpaceCellListImpl::get_structure_dimension(
    const Species& sp) const
{
    structure_container_type::const_iterator i(structures_.find(sp));
    if (i == structures_.end())
    {
        return Shape::THREE; // Default value (ex. for VACANT type)
    }
    return (*i).second->dimension();
}

} // ecell4
