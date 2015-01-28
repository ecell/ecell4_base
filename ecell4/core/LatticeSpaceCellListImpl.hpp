#ifndef __ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP
#define __ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP

#include "Context.hpp"
#include "MolecularType.hpp"
#include "VacantType.hpp"
// #include <cmath>
#include <sstream>

#include "LatticeSpace.hpp"

namespace ecell4
{

class LatticeSpaceCellListImpl
    : public LatticeSpaceBase
{
public:

    typedef LatticeSpaceBase base_type;

    typedef base_type::particle_info_type particle_info_type;
    typedef base_type::private_coordinate_type private_coordinate_type;
    typedef base_type::private_coordinate_type coordinate_type;

    typedef std::map<Species, MolecularType> spmap;
    // typedef std::vector<MolecularTypeBase*> voxel_container;

public:

    LatticeSpaceCellListImpl(
        const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic = true)
        : base_type(edge_lengths, voxel_radius), is_periodic_(is_periodic)
    {
        vacant_ = &(VacantType::getInstance());
        std::stringstream ss;
        ss << voxel_radius_;
        border_ = new MolecularType(Species("Border", ss.str(), "0"));
        periodic_ = new MolecularType(Species("Periodic", ss.str(), "0"));
    }

    virtual ~LatticeSpaceCellListImpl()
    {
        delete border_;
        delete periodic_;
    }

    virtual std::vector<Species> list_species() const
    {
        std::vector<Species> keys;
        keys.reserve(spmap_.size());
        for (spmap::const_iterator i(spmap_.begin());
            i != spmap_.end(); ++i)
        {
            keys.push_back((*i).first);
        }
        return keys;
    }

    virtual Integer num_voxels_exact(const Species& sp) const
    {
        spmap::const_iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            return 0;
        }
        const MolecularTypeBase* mt(&((*itr).second));
        return mt->size();
    }

    virtual Integer num_voxels(const Species& sp) const
    {
        Integer count(0);
        SpeciesExpressionMatcher sexp(sp);
        for (spmap::const_iterator itr(spmap_.begin());
                itr != spmap_.end(); ++itr)
        {
            if (sexp.match((*itr).first))
            {
                const MolecularTypeBase* mt(&((*itr).second));
                count += mt->size();
            }
        }
        return count;
    }

    virtual Integer num_voxels() const
    {
        Integer retval(0);
        for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
        {
            retval += (*itr).second.size();
        }
        return retval;
    }

    virtual bool has_voxel(const ParticleID& pid) const
    {
        for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
        {
            const MolecularTypeBase& mt((*itr).second);
            if (mt.is_vacant())
            {
                return false;
            }
            for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
            {
                if ((*vitr).second == pid)
                {
                    return true;
                }
            }
        }
        return false;
    }

    virtual std::vector<std::pair<ParticleID, Voxel> > list_voxels() const
    {
        std::vector<std::pair<ParticleID, Voxel> > retval;

        for (spmap::const_iterator itr(spmap_.begin()); itr != spmap_.end(); ++itr)
        {
            const MolecularTypeBase* mt(&((*itr).second));
            const std::string loc((mt->location()->is_vacant())
                ? "" : mt->location()->species().serial());
            const Species& sp(mt->species());
            for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
                itr != mt->end(); ++itr)
            {
                retval.push_back(std::make_pair(
                    (*itr).second,
                    Voxel(sp, private2coord((*itr).first), mt->radius(), mt->D(), loc)));
            }
        }
        return retval;
    }

    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const
    {
        SpeciesExpressionMatcher sexp(sp);
        std::vector<std::pair<ParticleID, Voxel> > retval;
        for (spmap::const_iterator itr(spmap_.begin());
                itr != spmap_.end(); ++itr)
        {
            if (!sexp.match((*itr).first))
            {
                continue;
            }

            const MolecularTypeBase* mt(&((*itr).second));
            const std::string loc((mt->location()->is_vacant())
                ? "" : mt->location()->species().serial());
            for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
                itr != mt->end(); ++itr)
            {
                retval.push_back(std::make_pair(
                    (*itr).second,
                    Voxel(sp, private2coord((*itr).first), mt->radius(), mt->D(), loc)));
            }
        }
        return retval;
    }

    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, Voxel> > retval;
        spmap::const_iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            return retval;
        }

        const MolecularTypeBase* mt(&((*itr).second));
        const std::string loc((mt->location()->is_vacant())
            ? "" : mt->location()->species().serial());
        for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
            itr != mt->end(); ++itr)
        {
            retval.push_back(std::make_pair(
                (*itr).second,
                Voxel(sp, private2coord((*itr).first), mt->radius(), mt->D(), loc)));
        }
        return retval;
    }

    virtual bool update_voxel_private(const Voxel& v)
    {
        return false;
    }

    virtual bool update_voxel_private(const ParticleID& pid, const Voxel& v)
    {
        return false;
    }

    virtual std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const
    {
        for (spmap::const_iterator i(spmap_.begin()); i != spmap_.end(); ++i)
        {
            const MolecularType& mt((*i).second);
            MolecularType::container_type::const_iterator j(mt.find(pid));
            if (j != mt.end())
            {
                const coordinate_type coord(private2coord((*j).first));
                const std::string loc((mt.location()->is_vacant())
                    ? "" : mt.location()->species().serial());
                return std::make_pair(
                    pid, Voxel((*i).first, coord, mt.radius(), mt.D(), loc));
            }
        }

        throw NotFound("voxel not found.");
    }

    virtual bool remove_voxel(const ParticleID& pid)
    {
        return true;
    }

    virtual bool remove_voxel_private(const private_coordinate_type& coord)
    {
        return true;
    }

    virtual bool move(const coordinate_type& from, const coordinate_type& to)
    {
        return true;
    }

    virtual const Particle particle_at(const coordinate_type& coord) const
    {
        return Particle(Species(""), coordinate2position(coord), voxel_radius_, 0); //XXX: DUMMY
    }

    virtual MolecularTypeBase* find_molecular_type(const Species& sp)
    {
        LatticeSpaceVectorImpl::spmap::iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            throw NotFound("MolecularType not found.");
        }
        return &((*itr).second);
    }

    virtual MolecularTypeBase* get_molecular_type(
        const private_coordinate_type& coord)
    {
        return vacant_;
    }

    virtual bool on_structure(const Voxel& v)
    {
        return true;
    }

    virtual std::pair<private_coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand)
    {
        const private_coordinate_type private_from(info.first);
        private_coordinate_type private_to(get_neighbor(private_from, nrand));
        return std::make_pair(private_to, false);
    }

    virtual Integer num_molecules() const
    {
        return LatticeSpace::num_molecules();
    }

    virtual Integer num_molecules(const Species& sp) const;

protected:

    std::pair<spmap::iterator, bool> __get_molecular_type(const Voxel& v);

protected:

    bool is_periodic_;

    spmap spmap_;
    // voxel_container voxels_;

    MolecularTypeBase* vacant_;
    MolecularTypeBase* border_;
    MolecularTypeBase* periodic_;
};

} // ecell4

#endif /* __ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP */
