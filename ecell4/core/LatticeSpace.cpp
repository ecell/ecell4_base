#include "Context.hpp"
#include "MolecularType.hpp"
#include "VacantType.hpp"
#include "LatticeSpace.hpp"
#include <cmath>
#include <sstream>

#ifdef WIN32_MSC
#include <boost/numeric/interval/detail/msvc_rounding_control.hpp>
#endif


namespace ecell4
{

#ifdef WIN32_MSC
double rint(const double x)
{
    return boost::numeric::interval_lib::detail::rint(x);
}

double round(const double x)
{
    return floor(x + 0.5);
}
#endif

LatticeSpaceVectorImpl::LatticeSpaceVectorImpl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const bool is_periodic) :
    base_type(edge_lengths, voxel_radius), is_periodic_(is_periodic)
{
    vacant_ = &(VacantType::getInstance());
    std::stringstream ss;
    ss << voxel_radius_;
    border_ = new MolecularType(Species("Border", ss.str(), "0"));
    periodic_ = new MolecularType(Species("Periodic", ss.str(), "0"));

    set_lattice_properties(is_periodic_);
}

LatticeSpaceVectorImpl::~LatticeSpaceVectorImpl()
{
    delete border_;
    delete periodic_;
}

/*
 * Coordinate transformations
 */

LatticeSpaceVectorImpl::coordinate_type LatticeSpaceVectorImpl::global2coord_(const Integer3& global,
        Integer col_size, Integer row_size, Integer layer_size) const
{
    return global.row +
        row_size * (global.col) +
        row_size * col_size * (global.layer);
}

const Integer3 LatticeSpaceVectorImpl::coord2global_(coordinate_type coord,
        Integer col_size, Integer row_size, Integer layer_size) const
{
    // Integer3 retval;
    // retval.col = coord / (row_size * layer_size);
    // retval.layer = (coord % (row_size * layer_size)) / row_size;
    // retval.row = (coord % (row_size * layer_size)) % row_size;
    // return retval;

    Integer3 retval;
    const Integer NUM_COLROW(row_size * col_size);
    const Integer LAYER(coord / NUM_COLROW);
    const Integer SURPLUS(coord - LAYER * NUM_COLROW);
    const Integer COL(SURPLUS / row_size);
    retval.col = COL;
    retval.layer = LAYER;
    retval.row = SURPLUS - COL * row_size;
    return retval;
}

// LatticeSpaceVectorImpl::coordinate_type LatticeSpaceVectorImpl::global2coord(const Integer3& global) const
// {
//     return global2coord_(global, col_size(), row_size(), layer_size());
// }

// LatticeSpaceVectorImpl::private_coordinate_type LatticeSpaceVectorImpl::global2private_coord(
//         const Integer3& global) const
// {
//     Integer3 modified_global;
//     modified_global.col = global.col + 1;
//     modified_global.row = global.row + 1;
//     modified_global.layer = global.layer + 1;
//     return global2coord_(modified_global, col_size_, row_size_, layer_size_);
// }

// Integer3 LatticeSpaceVectorImpl::coord2global(coordinate_type coord) const
// {
//     return coord2global_(coord, col_size(), row_size(), layer_size());
// }

// Integer3 LatticeSpaceVectorImpl::private_coord2global(
//     const private_coordinate_type& private_coord) const
// {
//     Integer3 retval(private_coord2private_global(private_coord));
//     retval.col--;
//     retval.row--;
//     retval.layer--;
//     return retval;
// }

// Real3 LatticeSpaceVectorImpl::global2position(const Integer3& global) const
// {
//     //the center point of a voxel
//     Real3 position;
//     position[0] = global.col * HCP_X;
//     position[1] = (global.col % 2) * HCP_L + HCP_Y * global.layer;
//     position[2] = (global.row * 2 + (global.layer + global.col) % 2)
//         * voxel_radius_;
//     return position;
// }

// Integer3 LatticeSpaceVectorImpl::position2global(const Real3& pos) const
// {
//     Integer3 global;
//     global.col = round(pos[0] / HCP_X);
//     global.layer = round((pos[1] - (global.col % 2) * HCP_L) / HCP_Y);
//     global.row = round((pos[2] / voxel_radius_
//                 - ((global.layer + global.col) % 2)) / 2);
//     return global;
// }

// Real3 LatticeSpaceVectorImpl::coordinate2position(const coordinate_type& coord) const
// {
//     return private2position(coord2private(coord));
// }
// 
// LatticeSpaceVectorImpl::coordinate_type LatticeSpaceVectorImpl::position2coordinate(
//         const Real3& pos) const
// {
//     return global2coord(position2global(pos));
// }
// 
// Real3 LatticeSpaceVectorImpl::private2position(
//     const private_coordinate_type& private_coord) const
// {
//     return global2position(private_coord2global(private_coord));
// }
// 
// LatticeSpaceVectorImpl::private_coordinate_type LatticeSpaceVectorImpl::position2private(
//         const Real3& pos) const
// {
//     return global2private_coord(position2global(pos));
// }
// 
// 
// LatticeSpaceVectorImpl::private_coordinate_type LatticeSpaceVectorImpl::coord2private(
//     const coordinate_type& coord) const
// {
//     return global2private_coord(coord2global(coord));
// }
// 
// LatticeSpaceVectorImpl::coordinate_type LatticeSpaceVectorImpl::private2coord(
//     const private_coordinate_type& private_coord) const
// {
//     return global2coord(private_coord2global(private_coord));
// }

LatticeSpaceVectorImpl::private_coordinate_type LatticeSpaceVectorImpl::get_neighbor(
        private_coordinate_type private_coord, Integer nrand) const
{
    const Integer NUM_COLROW(col_size_ * row_size_);
    const Integer NUM_ROW(row_size_);
    const bool odd_col(((private_coord % NUM_COLROW) / NUM_ROW) & 1);
    const bool odd_lay((private_coord / NUM_COLROW) & 1);

    switch(nrand)
    {
    case 1:
        return private_coord + 1;
    case 2:
        return private_coord + (odd_col ^ odd_lay) - NUM_ROW - 1;
    case 3:
        return private_coord + (odd_col ^ odd_lay) - NUM_ROW;
    case 4:
        return private_coord + (odd_col ^ odd_lay) + NUM_ROW - 1;
    case 5:
        return private_coord + (odd_col ^ odd_lay) + NUM_ROW;
    case 6:
        return private_coord - (2 * odd_col - 1) * NUM_COLROW - NUM_ROW;
    case 7:
        return private_coord - (2 * odd_col - 1) * NUM_COLROW + NUM_ROW;
    case 8:
        return private_coord + (odd_col ^ odd_lay) - NUM_COLROW - 1;
    case 9:
        return private_coord + (odd_col ^ odd_lay) - NUM_COLROW;
    case 10:
        return private_coord + (odd_col ^ odd_lay) + NUM_COLROW - 1;
    case 11:
        return private_coord + (odd_col ^ odd_lay) + NUM_COLROW;
    }
    return private_coord - 1;
}

void LatticeSpaceVectorImpl::set_lattice_properties(const bool is_periodic)
{
    // HCP_L = voxel_radius_ / sqrt(3.0);
    // HCP_X = voxel_radius_ * sqrt(8.0 / 3.0); //Lx
    // HCP_Y = voxel_radius_ * sqrt(3.0); //Ly

    // const Real lengthX = edge_lengths_[0];
    // const Real lengthY = edge_lengths_[1];
    // const Real lengthZ = edge_lengths_[2];

    // row_size_ = (Integer)rint((lengthZ/2)/voxel_radius_) + 2;
    // layer_size_ = (Integer)rint(lengthY/HCP_Y) + 2;
    // col_size_ = (Integer)rint(lengthX/HCP_X) + 2;

    private_coordinate_type voxel_size(
            col_size_ * row_size_ * layer_size_);
    voxels_.reserve(voxel_size);
    for (private_coordinate_type coord(0); coord < voxel_size; ++coord)
    {
        if (!is_inside(coord))
        {
            if (is_periodic)
                voxels_.push_back(periodic_);
            else
                voxels_.push_back(border_);
        }
        else
        {
            voxels_.push_back(vacant_);
        }
    }
}

Integer LatticeSpaceVectorImpl::num_species() const
{
    return spmap_.size();
}

Integer LatticeSpaceVectorImpl::num_molecules(const Species& sp) const
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

bool LatticeSpaceVectorImpl::has_voxel(const ParticleID& pid) const
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

// std::vector<std::pair<ParticleID, Particle> >
//     LatticeSpaceVectorImpl::list_particles() const
// {
//     std::vector<std::pair<ParticleID, Particle> > retval;
//     for (spmap::const_iterator itr(spmap_.begin());
//             itr != spmap_.end(); ++itr)
//     {
//         const MolecularTypeBase& mt((*itr).second);
//         if (mt.is_vacant())
//         {
//             continue;
//         }
// 
//         for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
//                 vitr != mt.end(); ++vitr)
//         {
//             retval.push_back(std::pair<ParticleID, Particle>(
//                         (*vitr).second, particle_at((*vitr).first)));
//         }
//     }
// 
//     return retval;
// }
// 
// std::vector<std::pair<ParticleID, Particle> >
//     LatticeSpaceVectorImpl::list_particles_exact(const Species& sp) const
// {
//     std::vector<std::pair<ParticleID, Particle> > retval;
//     spmap::const_iterator itr(spmap_.find(sp));
//     if (itr != spmap_.end())
//     {
//         const MolecularTypeBase& mt((*itr).second);
//         if (mt.is_vacant())
//         {
//             return retval;
//         }
// 
//         for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
//                 vitr != mt.end(); ++vitr)
//         {
//             retval.push_back(std::pair<ParticleID, Particle>(
//                         (*vitr).second, particle_at((*vitr).first)));
//         }
//     }
//     return retval;
// }
// 
// std::vector<std::pair<ParticleID, Particle> >
//     LatticeSpaceVectorImpl::list_particles(const Species& sp) const
// {
//     SpeciesExpressionMatcher sexp(sp);
//     std::vector<std::pair<ParticleID, Particle> > retval;
//     for (spmap::const_iterator itr(spmap_.begin());
//             itr != spmap_.end(); ++itr)
//     {
//         if (!sexp.match((*itr).first))
//         {
//             continue;
//         }
// 
//         const MolecularTypeBase& mt((*itr).second);
//         if (mt.is_vacant())
//         {
//             return retval;
//         }
// 
//         for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
//                 vitr != mt.end(); ++vitr)
//         {
//             retval.push_back(std::pair<ParticleID, Particle>(
//                         (*vitr).second, particle_at((*vitr).first)));
//         }
//     }
//     return retval;
// }

std::pair<ParticleID, Voxel>
LatticeSpaceVectorImpl::get_voxel(const ParticleID& pid) const
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
            return std::make_pair(pid, Voxel((*i).first, coord, mt.radius(), mt.D(), loc));
        }
    }

    throw NotFound("voxel not found.");
}

bool LatticeSpaceVectorImpl::update_structure(const Particle& p)
{
    //XXX: Particle does not have a location.
    Voxel v(p.species(), position2private(p.position()), p.radius(), p.D());
    return update_voxel_private(ParticleID(), v);
}

/*
 * original methods
 */

std::vector<Species> LatticeSpaceVectorImpl::list_species() const
{
    std::vector<Species> keys;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        keys.push_back((*itr).first);
    }
    return keys;
}

const Species& LatticeSpaceVectorImpl::find_species(std::string name) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        if ((*itr).first.serial() == name)
        {
            return (*itr).first;
        }
    }
    throw NotFound(name);
}

std::vector<LatticeSpaceVectorImpl::coordinate_type>
    LatticeSpaceVectorImpl::list_coords_exact(const Species& sp) const
{
    std::vector<coordinate_type> retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        return retval;
    }

    const MolecularTypeBase* mt(&((*itr).second));

    for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
            itr != mt->end(); ++itr)
    {
        retval.push_back(private2coord((*itr).first));
    }
    return retval;
}

std::vector<LatticeSpaceVectorImpl::coordinate_type> LatticeSpaceVectorImpl::list_coords(const Species& sp) const
{
    std::vector<coordinate_type> retval;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        if (!spmatch(sp, (*itr).first))
        {
            continue;
        }

        const MolecularTypeBase* mt(&((*itr).second));

        for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
                itr != mt->end(); ++itr)
        {
            retval.push_back(private2coord((*itr).first));
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Voxel> >
LatticeSpaceVectorImpl::list_voxels() const
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

std::vector<std::pair<ParticleID, Voxel> >
LatticeSpaceVectorImpl::list_voxels_exact(const Species& sp) const
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

std::vector<std::pair<ParticleID, Voxel> >
LatticeSpaceVectorImpl::list_voxels(const Species& sp) const
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

std::pair<LatticeSpaceVectorImpl::spmap::iterator, bool>
LatticeSpaceVectorImpl::__get_molecular_type(const Voxel& v)
{
    LatticeSpaceVectorImpl::spmap::iterator itr(spmap_.find(v.species()));
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
            std::pair<LatticeSpaceVectorImpl::spmap::iterator, bool> locval(
                spmap_.insert(LatticeSpaceVectorImpl::spmap::value_type(locsp, locmt)));
            if (!locval.second)
            {
                throw AlreadyExists(
                    "never reach here. find_molecular_type seems wrong.");
            }

            location = &((*locval.first).second);
        }
    }

    MolecularType mt(v.species(), location, v.radius(), v.D());
    std::pair<LatticeSpaceVectorImpl::spmap::iterator, bool> retval(
        spmap_.insert(LatticeSpaceVectorImpl::spmap::value_type(v.species(), mt)));
    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval;
}

MolecularTypeBase* LatticeSpaceVectorImpl::find_molecular_type(const Species& sp)
{
    LatticeSpaceVectorImpl::spmap::iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        throw NotFound("MolecularType not found.");
    }
    return &((*itr).second);
}

// MolecularTypeBase* LatticeSpaceVectorImpl::find_molecular_type(const std::string name)
// {
//     for (spmap::iterator itr(spmap_.begin());
//             itr != spmap_.end(); ++itr)
//     {
//         if ((*itr).first.serial() == name) {
//             return &((*itr).second);
//         }
//     }
// 
//     throw NotFound("MolecularType not found.");
// }

MolecularTypeBase* LatticeSpaceVectorImpl::get_molecular_type(const Voxel& v)
{
    return &((*(__get_molecular_type(v).first)).second);
}

bool LatticeSpaceVectorImpl::on_structure(const Voxel& v)
{
    // return get_molecular_type(v.coordinate()) != get_molecular_type(v)->location();
    return voxels_.at(v.coordinate()) != get_molecular_type(v)->location();
}

// bool LatticeSpaceVectorImpl::register_species(const Species& sp)
// {
//     return __get_molecular_type(sp).second;
// }

/*
 * Protected functions
 */

LatticeSpaceVectorImpl::private_coordinate_type LatticeSpaceVectorImpl::get_coord(const ParticleID& pid) const
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
                return (*vitr).first;
            }
        }
    }
    return -1;
}

MolecularTypeBase* LatticeSpaceVectorImpl::get_molecular_type(const private_coordinate_type& coord)
{
    // return voxels_.at(coord2private(coord));
    return voxels_.at(coord);
}

// bool LatticeSpaceVectorImpl::has_species_exact(const Species& sp) const
// {
//     return spmap_.find(sp) != spmap_.end();
// }

bool LatticeSpaceVectorImpl::has_species(const Species& sp) const
{
    return spmap_.find(sp) != spmap_.end();
}

bool LatticeSpaceVectorImpl::remove_voxel(const ParticleID& pid)
{
    for (spmap::iterator i(spmap_.begin()); i != spmap_.end(); ++i)
    {
        MolecularType& mt((*i).second);
        MolecularType::container_type::const_iterator j(mt.find(pid));
        if (j != mt.end())
        {
            const private_coordinate_type coord((*j).first);
            if (!mt.removeVoxel(coord))
            {
                return false;
            }

            voxel_container::iterator itr(voxels_.begin() + coord);
            (*itr) = mt.location();
            // mt.location()->addVoxel(particle_info_type(coord, ParticleID()));
            mt.location()->add_voxel_without_checking(
                particle_info_type(coord, ParticleID()));
            return true;
        }
    }
    return false;
}

bool LatticeSpaceVectorImpl::remove_voxel_private(const private_coordinate_type& coord)
{
    voxel_container::iterator itr(voxels_.begin() + coord);
    MolecularTypeBase* mt(*itr);
    if (mt->is_vacant())
    {
        return false;
    }
    if (mt->removeVoxel(coord))
    {
        (*itr) = mt->location();
        return true;
    }
    return false;
}

bool LatticeSpaceVectorImpl::move(const coordinate_type& from, const coordinate_type& to)
{
    const coordinate_type private_from(coord2private(from));
    const coordinate_type private_to(coord2private(to));
    return move_(private_from, private_to).second;
}

std::pair<LatticeSpaceVectorImpl::coordinate_type, bool> LatticeSpaceVectorImpl::move_to_neighbor(
        private_coordinate_type coord, Integer nrand)
{
    const private_coordinate_type neighbor(get_neighbor(coord, nrand));
    return move_(coord, neighbor);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool> LatticeSpaceVectorImpl::move_to_neighbor(
    particle_info_type& info, Integer nrand)
{
    const coordinate_type neighbor(get_neighbor(info.first, nrand));
    return move_(info, neighbor);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool> LatticeSpaceVectorImpl::move_(
    private_coordinate_type private_from, private_coordinate_type private_to)
{
    if (private_from == private_to)
    {
        // std::cerr << " from == to ";
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }

    MolecularTypeBase* from_mt(voxels_.at(private_from));
    if (from_mt->is_vacant())
    {
        return std::pair<private_coordinate_type, bool>(private_from, true);
    }

    MolecularTypeBase* to_mt(voxels_.at(private_to));

    if (to_mt == border_)
    {
        // std::cerr << " to_mt is border ";
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }
    else if (to_mt == periodic_)
    {
        private_to = apply_boundary_(private_to);
        to_mt = voxels_.at(private_to);
    }

    if (to_mt != from_mt->location())
    {
        // std::cerr << " to_mt is " << to_mt->species().serial() << " ";
        return std::pair<private_coordinate_type, bool>(private_to, false);
    }

    MolecularTypeBase::container_type::iterator itr(from_mt->find(private_from));
    (*itr).first = private_to;
    voxel_container::iterator from_itr(voxels_.begin() + private_from);
    (*from_itr) = to_mt;
    // to_mt->removeVoxel(private_to);
    // to_mt->addVoxel(particle_info_type(private_from, ParticleID()));
    to_mt->replace_voxel(private_to, particle_info_type(private_from, ParticleID()));
    voxel_container::iterator to_itr(voxels_.begin() + private_to);
    (*to_itr) = from_mt;

    return std::pair<private_coordinate_type, bool>(private_to, true);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool> LatticeSpaceVectorImpl::move_(
    particle_info_type& info, private_coordinate_type private_to)
{
    const private_coordinate_type private_from(info.first);
    if (private_from == private_to)
    {
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }

    MolecularTypeBase* from_mt(voxels_.at(private_from));
    if (from_mt->is_vacant())
    {
        return std::pair<private_coordinate_type, bool>(private_from, true);
    }

    MolecularTypeBase* to_mt(voxels_.at(private_to));

    if (to_mt == border_)
    {
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }
    else if (to_mt == periodic_)
    {
        private_to = apply_boundary_(private_to);
        to_mt = voxels_.at(private_to);
    }

    if (to_mt != from_mt->location())
    {
        return std::pair<private_coordinate_type, bool>(private_to, false);
    }

    voxel_container::iterator from_itr(voxels_.begin() + private_from);
    (*from_itr) = to_mt;
    // to_mt->removeVoxel(private_to);
    // to_mt->addVoxel(particle_info_type(private_from, ParticleID()));
    to_mt->replace_voxel(private_to, particle_info_type(private_from, ParticleID()));

    info.first = private_to;

    voxel_container::iterator to_itr(voxels_.begin() + private_to);
    (*to_itr) = from_mt;

    return std::pair<private_coordinate_type, bool>(private_to, true);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool> LatticeSpaceVectorImpl::move_to_neighbor(
    MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
    particle_info_type& info, const Integer nrand)
{
    const private_coordinate_type private_from(info.first);
    private_coordinate_type private_to(get_neighbor(private_from, nrand));

    //XXX: assert(private_from != private_to);
    //XXX: assert(from_mt == voxels_[private_from]);
    //XXX: assert(from_mt != vacant_);

    MolecularTypeBase* to_mt(voxels_[private_to]);

    if (to_mt != loc) // (to_mt != from_mt->location())
    {
        if (to_mt == border_)
        {
            return std::make_pair(private_from, false);
        }
        else if (to_mt != periodic_)
        {
            return std::make_pair(private_to, false);
        }

        // to_mt == periodic_
        private_to = apply_boundary_(private_to);
        to_mt = voxels_[private_to];

        if (to_mt != loc)
        {
            return std::make_pair(private_to, false);
        }
    }

    voxels_[private_from] = to_mt;
    voxels_[private_to] = from_mt;
    info.first = private_to; //XXX: updating data

    if (to_mt != vacant_) // (!to_mt->is_vacant())
    {
        to_mt->replace_voxel(
            private_to, particle_info_type(private_from, ParticleID()));
    }
    return std::make_pair(private_to, true);
}

// std::vector<LatticeSpaceVectorImpl::private_coordinate_type>
// LatticeSpaceVectorImpl::get_neighbors(LatticeSpaceVectorImpl::private_coordinate_type coord) const
// {
//     std::vector<LatticeSpaceVectorImpl::private_coordinate_type> retval;
//     for (Integer i(0); i < 12; i++)
//     {
//         retval.push_back(get_neighbor(coord, i));
//     }
//     return retval;
// }

const Particle LatticeSpaceVectorImpl::particle_at_private(private_coordinate_type coord) const
{
    const MolecularTypeBase* ptr_mt(voxels_.at(coord));
    // const Species& sp = ptr_mt->species();
    // const Real3& pos = coordinate2position(coord);
    // const Real& radius = 0;
    // const Real& D = 0;
    // Particle particle(sp, pos, radius, D);
    // return particle;
    return Particle(ptr_mt->species(), coordinate2position(
                private2coord(coord)), ptr_mt->radius(), ptr_mt->D());
}

bool LatticeSpaceVectorImpl::is_in_range(coordinate_type coord) const
{
    return coord >= 0 && coord < row_size() * layer_size() * col_size();
}

bool LatticeSpaceVectorImpl::is_in_range_private(private_coordinate_type coord) const
{
    return coord >= 0 && coord < row_size_ * col_size_ * layer_size_;
}

bool LatticeSpaceVectorImpl::is_inside(private_coordinate_type coord) const
{
    const Integer3 global(private_coord2private_global(coord));
    return global.col > 0 && global.col < col_size_-1
        && global.row > 0 && global.row < row_size_-1
        && global.layer > 0 && global.layer < layer_size_-1;
}

/*
 * Coordinate transformations
 */

const Integer3 LatticeSpaceVectorImpl::private_coord2private_global(
        const private_coordinate_type private_coord) const
{
    return coord2global_(private_coord, col_size_, row_size_, layer_size_);
}

LatticeSpaceVectorImpl::private_coordinate_type LatticeSpaceVectorImpl::apply_boundary_(
        const private_coordinate_type& private_coord) const
{
    Integer3 global(private_coord2private_global(private_coord));

    global.col = (global.col - 1) % col_size();
    global.row = (global.row - 1) % row_size();
    global.layer = (global.layer - 1) % layer_size();

    global.col = global.col < 0 ? global.col + col_size() : global.col;
    global.row = global.row < 0 ? global.row + row_size() : global.row;
    global.layer = global.layer < 0 ? global.layer + layer_size() : global.layer;

    return global2private_coord(global);
}

Integer LatticeSpaceVectorImpl::num_voxels_exact(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        return 0;
    }
    const MolecularTypeBase* mt(&((*itr).second));
    return mt->size();
}

Integer LatticeSpaceVectorImpl::num_voxels(const Species& sp) const
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

Integer LatticeSpaceVectorImpl::num_voxels() const
{
    Integer retval(0);
    for (spmap::const_iterator itr(spmap_.begin()); itr != spmap_.end(); ++itr)
    {
        retval += (*itr).second.size();
    }
    return retval;
}

bool LatticeSpaceVectorImpl::update_voxel_private(const Voxel& v)
{
    const LatticeSpaceVectorImpl::private_coordinate_type coord(v.coordinate());
    MolecularTypeBase* src_mt(voxels_.at(coord));
    if (src_mt->is_vacant())
    {
        return false;
    }

    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    if (new_mt->is_vacant())
    {
        return false; // Vacant has no ParticleID. Call remove_voxel.
    }

    // MolecularTypeBase::particle_info_type info(*(src_mt->find(coord)));
    // src_mt->removeVoxel(coord);
    // new_mt->addVoxel(info);
    MolecularTypeBase::iterator position(src_mt->find(coord));
    new_mt->add_voxel_without_checking(*position);
    src_mt->remove_voxel(position);
    voxel_container::iterator itr(voxels_.begin() + coord);
    (*itr) = new_mt;
    return true;
}

bool LatticeSpaceVectorImpl::update_voxel_private(const ParticleID& pid, const Voxel& v)
{
    const LatticeSpaceVectorImpl::private_coordinate_type& to_coord(v.coordinate());
    if (!is_in_range_private(to_coord))
    {
        return false;
    }

    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    if (new_mt->is_vacant())
    {
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
        const LatticeSpaceVectorImpl::private_coordinate_type from_coord(get_coord(pid));
        if (from_coord != -1)
        {
            MolecularTypeBase* src_mt(voxels_.at(from_coord));
            src_mt->removeVoxel(from_coord);
            voxel_container::iterator itr(voxels_.begin() + from_coord);
            // voxels_.erase(itr);
            // voxels_.insert(itr, dest_mt);
            (*itr) = dest_mt;
            // dest_mt->addVoxel(particle_info_type(from_coord, ParticleID()));
            // dest_mt->add_voxel_without_checking(
            //     particle_info_type(from_coord, ParticleID()));
            // dest_mt->removeVoxel(to_coord);
            dest_mt->replace_voxel(to_coord, particle_info_type(from_coord, ParticleID()));
            new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
            itr = voxels_.begin() + to_coord;
            (*itr) = new_mt;
            return true;
        }
    }

    // new_mt->addVoxel(MolecularTypeBase::particle_info_type(to_coord, pid));
    new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
    voxel_container::iterator to_itr(voxels_.begin() + to_coord);
    (*to_itr) = new_mt;
    dest_mt->removeVoxel(to_coord);
    return true;
}

} // ecell4
