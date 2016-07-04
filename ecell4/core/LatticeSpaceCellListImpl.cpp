#include "LatticeSpaceCellListImpl.hpp"
#include "StructureType.hpp"


namespace ecell4
{

Integer LatticeSpaceCellListImpl::num_molecules(const Species& sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            if (!mt->with_voxels())
            {
                count += count_voxels(mt);
                // throw NotSupported(
                //     "num_molecules for MolecularType with no voxel"
                //     " is not supporeted now");
            }
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
    VoxelPool* dest_mt(get_molecular_type(to_coord));

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

    // const coordinate_type& to_coord(v.coordinate());
    // if (!is_in_range(to_coord))
    // {
    //     return false;
    // }

    // VoxelPool* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    // if (new_mt->is_vacant())
    // {
    //     // ???
    //     return false; // Vacant has no ParticleID.
    // }

    // VoxelPool* dest_mt(get_molecular_type(to_coord));
    // if (!dest_mt->is_vacant() && dest_mt != new_mt->location())
    // {
    //     if (dest_mt->species() == periodic_->species()
    //         || dest_mt->species() == border_->species())
    //     {
    //         throw NotSupported("The coordinate points a boundary.");
    //     }

    //     const ParticleID to_pid(dest_mt->find_particle_id(to_coord));
    //     if (pid == ParticleID() || to_pid != pid)
    //     {
    //         return false; // collision
    //     }

    //     // VoxelPool::container_type::const_iterator
    //     //     dest_itr(dest_mt->find(to_coord));
    //     // if (dest_itr == dest_mt->end())
    //     // {
    //     //     throw IllegalState(
    //     //         "MolecularTypaBase [" + dest_mt->species().serial()
    //     //         + "] doesn't contain a proper coordinate.");
    //     // }
    //     // else if (!(pid != ParticleID() && (*dest_itr).second == pid))
    //     // {
    //     //     return false; // collision
    //     // }
    // }

    // if (pid != ParticleID())
    // {
    //     const std::pair<VoxelPool*, coordinate_type>
    //         target(__get_coordinate(pid));
    //     const coordinate_type& from_coord(target.second);
    //     if (from_coord != -1)
    //     {
    //         VoxelPool* src_mt(target.first);
    //         src_mt->remove_voxel_if_exists(from_coord);
    //         dest_mt->replace_voxel(to_coord, coordinate_id_pair_type(ParticleID(), from_coord));
    //         new_mt->add_voxel_without_checking(coordinate_id_pair_type(pid, to_coord));

    //         if (!dest_mt->is_vacant())
    //         {
    //             update_matrix(from_coord, dest_mt);
    //             update_matrix(to_coord, new_mt);
    //         }
    //         else
    //         {
    //             update_matrix(from_coord, to_coord, new_mt);
    //         }
    //         return true;
    //     }
    // }

    // new_mt->add_voxel_without_checking(coordinate_id_pair_type(pid, to_coord));
    // dest_mt->remove_voxel_if_exists(to_coord);
    // update_matrix(to_coord, new_mt);
    // return true;
}

std::pair<LatticeSpaceCellListImpl::spmap::iterator, bool>
LatticeSpaceCellListImpl::__get_molecular_type(const Voxel& v)
{
    spmap::iterator itr(spmap_.find(v.species()));
    if (itr != spmap_.end())
    {
        return std::make_pair(itr, false);
    }

    VoxelPool* location;
    if (v.loc() == "")
    {
        location = vacant_;
    }
    else
    {
        const Species locsp(v.loc());
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
            boost::shared_ptr<MolecularType> locmt(
                new MolecularType(locsp, vacant_, voxel_radius_, 0));
            std::pair<spmap::iterator, bool> locval(
                spmap_.insert(spmap::value_type(locsp, locmt)));
            if (!locval.second)
            {
                throw AlreadyExists(
                    "never reach here. find_molecular_type seems wrong.");
            }

            location = (*locval.first).second.get();
        }
    }

    boost::shared_ptr<MolecularType> mt(
        new MolecularType(v.species(), location, v.radius(), v.D()));
    std::pair<spmap::iterator, bool> retval(
        spmap_.insert(spmap::value_type(v.species(), mt)));
    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval;
}

std::pair<VoxelPool*, LatticeSpaceCellListImpl::coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid)
{
    for (spmap::iterator itr(spmap_.begin());
        itr != spmap_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        for (MolecularType::const_iterator vitr(mt->begin());
            vitr != mt->end(); ++vitr)
        {
            if ((*vitr).pid == pid)
            {
                return std::pair<VoxelPool*, coordinate_type>(
                    mt.get(), (*vitr).coordinate);
            }
        }
    }
    return std::make_pair<VoxelPool*, coordinate_type>(
        NULL, -1); //XXX: a bit dirty way
}

std::pair<const VoxelPool*, LatticeSpaceCellListImpl::coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
        itr != spmap_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        for (VoxelPool::const_iterator vitr(mt->begin());
            vitr != mt->end(); ++vitr)
        {
            if ((*vitr).pid == pid)
            {
                return std::pair<const VoxelPool*, coordinate_type>(
                    mt.get(), (*vitr).coordinate);
            }
        }
    }
    return std::make_pair<const VoxelPool*, coordinate_type>(
        NULL, -1); //XXX: a bit dirty way
}

VoxelPool* LatticeSpaceCellListImpl::get_molecular_type(
    const LatticeSpaceCellListImpl::coordinate_type& coord)
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

    cell_type& cell(matrix_[coordinate2index(coord)]);
    if (cell.size() == 0)
    {
        return vacant_;
    }

    cell_type::iterator i(find_from_cell(coord, cell));
    if (i != cell.end())
    {
        return (*i).first;
    }

    return vacant_;
}

const VoxelPool* LatticeSpaceCellListImpl::get_molecular_type(
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

    // for (spmap::const_iterator itr(spmap_.begin());
    //     itr != spmap_.end(); ++itr)
    // {
    //     const VoxelPool& mt((*itr).second);
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
    return (*(__get_molecular_type(v).first)).second.get();
}

std::pair<LatticeSpaceCellListImpl::coordinate_type, bool>
    LatticeSpaceCellListImpl::move_to_neighbor(
        VoxelPool* const& from_mt, VoxelPool* const& loc,
        LatticeSpaceCellListImpl::coordinate_id_pair_type& info, const Integer nrand)
{
    const coordinate_type from(info.coordinate);
    coordinate_type to(get_neighbor(from, nrand));

    VoxelPool* to_mt(get_molecular_type(to));

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
        to_mt = get_molecular_type(to);

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

    // const coordinate_type from(info.first);
    // coordinate_type to(get_neighbor(from, nrand));
    // VoxelPool* to_mt(get_molecular_type(to));
    // if (to_mt != loc)
    // {
    //     if (to_mt == border_)
    //     {
    //         return std::make_pair(from, false);
    //     }
    //     else if (to_mt != periodic_)
    //     {
    //         return std::make_pair(to, false);
    //     }

    //     // to_mt == periodic_
    //     to = periodic_transpose(to);
    //     to_mt = get_molecular_type(to);
    //     if (to_mt != loc)
    //     {
    //         return std::make_pair(to, false);
    //     }
    // }

    // info.first = to;
    // if (to_mt != vacant_) // (!to_mt->is_vacant())
    // {
    //     to_mt->replace_voxel(
    //         to, coordinate_id_pair_type(ParticleID(), from));
    //     update_matrix(from, to_mt);
    //     update_matrix(to, from_mt);
    // }
    // else
    // {
    //     update_matrix(from, to, from_mt);
    // }
    // return std::make_pair(to, true);
}

bool LatticeSpaceCellListImpl::make_structure_type(
    const Species& sp, Shape::dimension_kind dimension, const std::string loc)
{
    spmap::iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
    {
        return false;
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
            std::pair<LatticeSpaceVectorImpl::spmap::iterator, bool> locval(
                spmap_.insert(LatticeSpaceVectorImpl::spmap::value_type(locsp, locmt)));
            if (!locval.second)
            {
                throw AlreadyExists(
                    "never reach here. find_molecular_type seems wrong.");
            }

            location = (*locval.first).second.get();
        }
    }

    boost::shared_ptr<MolecularType> mt(new StructureType(sp, location, voxel_radius_, dimension));
    std::pair<spmap::iterator, bool>
        retval(spmap_.insert(std::make_pair(sp, mt)));
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
