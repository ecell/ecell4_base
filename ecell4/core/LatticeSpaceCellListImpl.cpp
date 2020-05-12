#include "LatticeSpaceCellListImpl.hpp"
#include "StructureType.hpp"

namespace ecell4
{

Integer LatticeSpaceCellListImpl::num_molecules(const Species &sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    for (const auto &info : molecule_pools_)
    {
        const Integer cnt(sexp.count(info.first));
        if (cnt > 0)
        {
            const std::shared_ptr<MoleculePool> &vp(info.second);
            count += vp->size() * cnt;
        }
    }
    return count;
}

/*
 * Change the Species and coordinate of a Voxel with ParticleID, pid, to
 * species and coordinate respectively and return false.
 * If no Voxel with pid is found, create a new Voxel at
 * coordiante() and return true.
 */
bool LatticeSpaceCellListImpl::update_voxel(const ParticleID &pid,
                                            const Species &species,
                                            const coordinate_type to_coord)
{
    if (!is_in_range(to_coord))
    {
        throw NotSupported("Out of bounds");
    }

    std::shared_ptr<VoxelPool> new_vp(find_voxel_pool(species));
    std::shared_ptr<VoxelPool> dest_vp(get_voxel_pool_at(to_coord));

    if (dest_vp != new_vp->location())
    {
        throw NotSupported("Mismatch in the location.");
    }

    if (pid != ParticleID())
    {
        const std::pair<std::shared_ptr<VoxelPool>, coordinate_type> target(
            __get_coordinate(pid));
        const coordinate_type &from_coord(target.second);
        if (from_coord != -1)
        {
            // move
            target.first->remove_voxel_if_exists(from_coord);

            // XXX: use location?
            dest_vp->replace_voxel(to_coord, from_coord);

            new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));

            if (!dest_vp->is_vacant())
            {
                update_matrix(from_coord, dest_vp);
                update_matrix(to_coord, new_vp);
            }
            else
            {
                update_matrix(from_coord, to_coord, new_vp);
            }
            return true;
        }
    }

    // new
    dest_vp->remove_voxel_if_exists(to_coord);

    new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));
    update_matrix(to_coord, new_vp);
    return true;
}

bool LatticeSpaceCellListImpl::add_voxel(const Species &sp,
                                         const ParticleID &pid,
                                         const coordinate_type &coord)
{
    std::shared_ptr<VoxelPool> vpool(find_voxel_pool(sp));
    std::shared_ptr<VoxelPool> location(get_voxel_pool_at(coord));

    if (vpool->location() != location)
        return false;

    location->remove_voxel_if_exists(coord);
    vpool->add_voxel(coordinate_id_pair_type(pid, coord));
    update_matrix(coord, vpool);

    return true;
}

std::pair<std::shared_ptr<VoxelPool>,
          LatticeSpaceCellListImpl::coordinate_type>
LatticeSpaceCellListImpl::__get_coordinate(const ParticleID &pid)
{
    for (auto &info : molecule_pools_)
    {
        const std::shared_ptr<MoleculePool> &vp(info.second);
        MoleculePool::container_type::const_iterator j(vp->find(pid));
        if (j != vp->end())
        {
            return std::make_pair(vp, (*j).coordinate);
        }
    }
    return std::make_pair(std::shared_ptr<VoxelPool>(),
                          -1); // XXX: a bit dirty way
}

std::pair<std::shared_ptr<const VoxelPool>,
          LatticeSpaceCellListImpl::coordinate_type>
LatticeSpaceCellListImpl::__get_coordinate(const ParticleID &pid) const
{
    for (const auto &info : molecule_pools_)
    {
        const std::shared_ptr<MoleculePool> &vp(info.second);
        MoleculePool::container_type::const_iterator j(vp->find(pid));
        if (j != vp->end())
        {
        }
    }
    return std::make_pair(std::shared_ptr<VoxelPool>(),
                          -1); // XXX: a bit dirty way
}

std::shared_ptr<VoxelPool> LatticeSpaceCellListImpl::get_voxel_pool_at(
    const LatticeSpaceCellListImpl::coordinate_type &coord) const
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
    //     VoxelPool& vp((*itr).second);
    //     if (vp.is_vacant())
    //     {
    //         continue;
    //     }

    //     VoxelPool::container_type::const_iterator j(vp.find(coord));
    //     if (j != vp.end())
    //     {
    //         return (&vp);
    //     }
    // }

    const cell_type &cell(matrix_[coordinate2index(coord)]);
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

void LatticeSpaceCellListImpl::add_structure(
    const Species &sp, const std::shared_ptr<const Shape> &s,
    const std::string loc)
{
    make_structure_type(sp, loc);

    structure_container_type::const_iterator i(structures_.find(sp));
    if (i != structures_.end())
    {
        throw NotSupported("not supported yet.");
    }
    structures_.insert(std::make_pair(sp, s));
}

const std::shared_ptr<const Shape> &
LatticeSpaceCellListImpl::get_structure(const Species &sp) const
{
    structure_container_type::const_iterator i(structures_.find(sp));
    if (i == structures_.end())
    {
        throw NotFound("not found.");
    }
    return (*i).second;
}

const Shape::dimension_kind
LatticeSpaceCellListImpl::get_structure_dimension(const Species &sp) const
{
    structure_container_type::const_iterator i(structures_.find(sp));
    if (i == structures_.end())
    {
        return Shape::THREE; // Default value (ex. for VACANT type)
    }
    return (*i).second->dimension();
}

} // namespace ecell4
