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

} // ecell4
