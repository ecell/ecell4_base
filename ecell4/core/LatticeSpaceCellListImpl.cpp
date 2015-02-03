#include "LatticeSpaceCellListImpl.hpp"


namespace ecell4
{

Integer LatticeSpaceCellListImpl::num_molecules(const Species& sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase* mt(&((*itr).second));
        count += mt->size() * sexp.count((*itr).first);
    }
    return count;
}

bool LatticeSpaceCellListImpl::update_voxel_private(const ParticleID& pid, const Voxel& v)
{
    const private_coordinate_type& to_coord(v.coordinate());
    if (!is_in_range_private(to_coord))
    {
        return false;
    }

    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    if (new_mt->is_vacant())
    {
        // ???
        return false; // Vacant has no ParticleID.
    }

    MolecularTypeBase* dest_mt(get_molecular_type(to_coord));
    if (!dest_mt->is_vacant() && dest_mt != new_mt->location())
    {
        if (dest_mt->species() == periodic_->species()
            || dest_mt->species() == border_->species())
        {
            throw NotSupported("The coordinate points a boundary.");
        }

        MolecularTypeBase::container_type::const_iterator
            dest_itr(dest_mt->find(to_coord));
        if (dest_itr == dest_mt->end())
        {
            throw IllegalState(
                "MolecularTypaBase [" + dest_mt->species().serial()
                + "] doesn't contain a proper coordinate.");
        }
        else if (!(pid != ParticleID() && (*dest_itr).second == pid))
        {
            return false; // collision
        }
    }

    if (pid != ParticleID())
    {
        const std::pair<MolecularTypeBase*, private_coordinate_type>
            target(__get_coordinate(pid));
        const private_coordinate_type& from_coord(target.second);
        if (from_coord != -1)
        {
            MolecularTypeBase* src_mt(target.first);
            src_mt->removeVoxel(from_coord);
            dest_mt->replace_voxel(to_coord, particle_info_type(from_coord, ParticleID()));
            new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));

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

    new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
    dest_mt->removeVoxel(to_coord);
    update_matrix(to_coord, new_mt);
    return true;
}

std::pair<LatticeSpaceCellListImpl::spmap::iterator, bool>
LatticeSpaceCellListImpl::__get_molecular_type(const Voxel& v)
{
    spmap::iterator itr(spmap_.find(v.species()));
    if (itr != spmap_.end())
    {
        return std::make_pair(itr, false);
    }

    MolecularTypeBase* location;
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
            // XXX: A MolecularTypeBase for the structure (location) must be allocated
            // XXX: before the allocation of a Species on the structure.
            // XXX: The MolecularTypeBase cannot be automatically allocated at the time
            // XXX: because its MoleculeInfo is unknown.
            // XXX: LatticeSpaceVectorImpl::load will raise a problem about this issue.
            // XXX: In this implementation, the MolecularTypeBase for a structure is
            // XXX: created with default arguments.
            MolecularType locmt(locsp, vacant_, voxel_radius_, 0);
            std::pair<spmap::iterator, bool> locval(
                spmap_.insert(spmap::value_type(locsp, locmt)));
            if (!locval.second)
            {
                throw AlreadyExists(
                    "never reach here. find_molecular_type seems wrong.");
            }

            location = &((*locval.first).second);
        }
    }

    MolecularType mt(v.species(), location, v.radius(), v.D());
    std::pair<spmap::iterator, bool> retval(
        spmap_.insert(spmap::value_type(v.species(), mt)));
    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval;
}

std::pair<MolecularTypeBase*, LatticeSpaceCellListImpl::private_coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid)
{
    for (spmap::iterator itr(spmap_.begin());
        itr != spmap_.end(); ++itr)
    {
        MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            continue;
        }

        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
            vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return std::pair<MolecularTypeBase*, private_coordinate_type>(
                    &mt, (*vitr).first);
            }
        }
    }
    return std::make_pair<MolecularTypeBase*, private_coordinate_type>(NULL, -1);
}

std::pair<const MolecularTypeBase*, LatticeSpaceCellListImpl::private_coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
        itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            continue;
        }

        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
            vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return std::pair<const MolecularTypeBase*, private_coordinate_type>(
                    &mt, (*vitr).first);
            }
        }
    }
    return std::make_pair<const MolecularTypeBase*, private_coordinate_type>(NULL, -1);
}

} // ecell4
