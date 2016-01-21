#include "Context.hpp"
#include "MolecularType.hpp"
#include "VacantType.hpp"
#include "StructureType.hpp"
#include "LatticeSpace.hpp"
#include <cmath>
#include <sstream>
#include <algorithm>

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

bool LatticeSpace::can_move(const private_coordinate_type& src,
        const private_coordinate_type& dest) const
{
    return false;
}

bool LatticeSpace::make_structure_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    return false;
}

void LatticeSpaceBase::set_lattice_properties(const bool is_periodic)
{
    //XXX: derived from SpatiocyteStepper::setLatticeProperties()
    HCP_L = voxel_radius_ / sqrt(3.0);
    HCP_X = voxel_radius_ * sqrt(8.0 / 3.0); // Lx
    HCP_Y = voxel_radius_ * sqrt(3.0); // Ly

    const Real lengthX = edge_lengths_[0];
    const Real lengthY = edge_lengths_[1];
    const Real lengthZ = edge_lengths_[2];

    col_size_ = (Integer)rint(lengthX / HCP_X) + 1;
    layer_size_ = (Integer)rint(lengthY / HCP_Y) + 1;
    row_size_ = (Integer)rint((lengthZ / 2) / voxel_radius_) + 1;

    if (is_periodic)
    {
        // The number of voxels in each axis must be even for a periodic boundary.
        col_size_ = (col_size_ % 2 == 0 ? col_size_ : col_size_ + 1);
        layer_size_ = (layer_size_ % 2 == 0 ? layer_size_ : layer_size_ + 1);
        row_size_ = (row_size_ % 2 == 0 ? row_size_ : row_size_ + 1);
    }

    row_size_ += 2;
    layer_size_ += 2;
    col_size_ += 2;
}

Integer3 LatticeSpaceBase::position2global(const Real3& pos) const
{
    const Integer col(round(pos[0] / HCP_X));
    const Integer layer(round((pos[1] - (col % 2) * HCP_L) / HCP_Y));
    const Integer row(round(
        (pos[2] / voxel_radius_ - ((layer + col) % 2)) / 2));
    const Integer3 global(col, row, layer);
    return global;
}

LatticeSpaceVectorImpl::LatticeSpaceVectorImpl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const bool is_periodic) :
    base_type(edge_lengths, voxel_radius, is_periodic), is_periodic_(is_periodic)
{
    vacant_ = &(VacantType::getInstance());
    std::stringstream ss;
    ss << voxel_radius_;
    border_ = new MolecularType(Species("Border", ss.str(), "0"));
    periodic_ = new MolecularType(Species("Periodic", ss.str(), "0"));

    initialize_voxels(is_periodic_);
}

LatticeSpaceVectorImpl::~LatticeSpaceVectorImpl()
{
    delete border_;
    delete periodic_;
}

void LatticeSpaceVectorImpl::initialize_voxels(const bool is_periodic)
{
    const private_coordinate_type voxel_size(
        col_size_ * row_size_ * layer_size_);
    // std::cout << "voxel_size = " << voxel_size << std::endl;

    spmap_.clear();
    voxels_.clear();
    voxels_.reserve(voxel_size);
    for (private_coordinate_type coord(0); coord < voxel_size; ++coord)
    {
        if (!is_inside(coord))
        {
            if (is_periodic)
            {
                voxels_.push_back(periodic_);
            }
            else
            {
                voxels_.push_back(border_);
            }
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
        const boost::shared_ptr<MolecularType>& mt((*itr).second);

        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            if (!mt->with_voxels())
            {
                // throw NotSupported(
                //     "num_molecules for MolecularType with no voxel"
                //     " is not supporeted now");
                count += count_voxels(mt) * cnt;
            }
            else
            {
                count += mt->size() * cnt;
            }
        }
    }
    return count;
}

bool LatticeSpaceVectorImpl::has_voxel(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);

        if (mt->is_vacant())
        {
            return false;
        }
        for (MolecularType::const_iterator vitr(mt->begin());
            vitr != mt->end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return true;
            }
        }
    }
    return false;
}

std::pair<ParticleID, Voxel>
LatticeSpaceVectorImpl::get_voxel(const ParticleID& pid) const
{
    for (spmap::const_iterator i(spmap_.begin()); i != spmap_.end(); ++i)
    {
        const boost::shared_ptr<MolecularType>& mt((*i).second);

        MolecularType::container_type::const_iterator j(mt->find(pid));
        if (j != mt->end())
        {
            const coordinate_type coord(private2coord((*j).first));
            const std::string loc((mt->location()->is_vacant())
                ? "" : mt->location()->species().serial());
            return std::make_pair(pid,
                Voxel((*i).first, coord, mt->radius(), mt->D(), loc));
        }
    }

    throw NotFound("voxel not found.");
}

std::pair<ParticleID, Voxel>
LatticeSpaceVectorImpl::get_voxel(const coordinate_type& coord) const
{
    const private_coordinate_type private_coord(coord2private(coord));
    const MolecularTypeBase* mt(voxels_[private_coord]);
    const std::string loc((mt->location()->is_vacant())
        ? "" : mt->location()->species().serial());
    if (mt->with_voxels())
    {
        return std::make_pair(mt->find_particle_id(private_coord),
            Voxel(mt->species(), coord, mt->radius(), mt->D(), loc));
    }
    else
    {
        return std::make_pair(ParticleID(),
            Voxel(mt->species(), coord, mt->radius(), mt->D(), loc));
    }
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

    const boost::shared_ptr<MolecularType>& mt((*itr).second);

    for (MolecularType::const_iterator itr(mt->begin()); itr != mt->end(); ++itr)
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

        const boost::shared_ptr<MolecularType>& mt((*itr).second);

        for (MolecularType::const_iterator itr(mt->begin());
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
        const boost::shared_ptr<MolecularType>& mt((*itr).second);

        const std::string loc((mt->location()->is_vacant())
            ? "" : mt->location()->species().serial());
        const Species& sp(mt->species());
        if (mt->with_voxels())
        {
            for (MolecularType::const_iterator i(mt->begin());
                i != mt->end(); ++i)
            {
                retval.push_back(std::make_pair(
                    (*i).second,
                    Voxel(sp, private2coord((*i).first), mt->radius(), mt->D(), loc)));
            }
        }
        else
        {
            for (voxel_container::const_iterator i(voxels_.begin());
                i != voxels_.end(); ++i)
            {
                if (*i != mt.get())
                {
                    continue;
                }

                const private_coordinate_type
                    private_coord(std::distance(voxels_.begin(), i));
                retval.push_back(std::make_pair(
                    ParticleID(),
                    Voxel(sp, private2coord(private_coord), mt->radius(), mt->D(), loc)));
            }
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

    const boost::shared_ptr<MolecularType>& mt((*itr).second);
    const std::string loc((mt->location()->is_vacant())
        ? "" : mt->location()->species().serial());
    if (mt->with_voxels())
    {
        for (MolecularType::const_iterator i(mt->begin());
            i != mt->end(); ++i)
        {
            retval.push_back(std::make_pair(
                (*i).second,
                Voxel(sp, private2coord((*i).first), mt->radius(), mt->D(), loc)));
        }
    }
    else
    {
        for (voxel_container::const_iterator i(voxels_.begin());
            i != voxels_.end(); ++i)
        {
            if (*i != mt.get())
            {
                continue;
            }

            const private_coordinate_type
                private_coord(std::distance(voxels_.begin(), i));
            retval.push_back(std::make_pair(
                ParticleID(),
                Voxel(sp, private2coord(private_coord), mt->radius(), mt->D(), loc)));
        }
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

        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        const std::string loc((mt->location()->is_vacant())
            ? "" : mt->location()->species().serial());
        if (mt->with_voxels())
        {
            for (MolecularType::const_iterator i(mt->begin());
                i != mt->end(); ++i)
            {
                retval.push_back(std::make_pair(
                    (*i).second,
                    Voxel(sp, private2coord((*i).first), mt->radius(), mt->D(), loc)));
            }
        }
        else
        {
            for (voxel_container::const_iterator i(voxels_.begin());
                i != voxels_.end(); ++i)
            {
                if (*i != mt.get())
                {
                    continue;
                }

                const private_coordinate_type
                    private_coord(std::distance(voxels_.begin(), i));
                retval.push_back(std::make_pair(
                    ParticleID(),
                    Voxel(sp, private2coord(private_coord), mt->radius(), mt->D(), loc)));
            }
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

    boost::shared_ptr<MolecularType>
        mt(new MolecularType(v.species(), location, v.radius(), v.D()));
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
    return (*itr).second.get(); //XXX: Raw pointer was thrown.
}

const MolecularTypeBase* LatticeSpaceVectorImpl::find_molecular_type(const Species& sp) const
{
    LatticeSpaceVectorImpl::spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        throw NotFound("MolecularType not found.");
    }
    return (*itr).second.get(); //XXX: Raw pointer was thrown.
}

MolecularTypeBase* LatticeSpaceVectorImpl::get_molecular_type(const Voxel& v)
{
    return (*(__get_molecular_type(v).first)).second.get(); //XXX: Raw pointer was thrown.
}

bool LatticeSpaceVectorImpl::on_structure(const Voxel& v)
{
    // return get_molecular_type(v.coordinate()) != get_molecular_type(v)->location();
    return voxels_.at(v.coordinate()) != get_molecular_type(v)->location();
}

/*
 * Protected functions
 */

LatticeSpaceVectorImpl::private_coordinate_type LatticeSpaceVectorImpl::get_coord(
    const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        for (MolecularType::const_iterator vitr(mt->begin());
            vitr != mt->end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return (*vitr).first;
            }
        }
    }
    return -1; //XXX: a bit dirty way
}

MolecularTypeBase* LatticeSpaceVectorImpl::get_molecular_type(
    const private_coordinate_type& coord)
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
        const boost::shared_ptr<MolecularType>& mt((*i).second);
        MolecularType::const_iterator j(mt->find(pid));
        if (j != mt->end())
        {
            const private_coordinate_type coord((*j).first);
            if (!mt->remove_voxel_if_exists(coord))
            {
                return false;
            }

            voxel_container::iterator itr(voxels_.begin() + coord);
            (*itr) = mt->location();
            mt->location()->add_voxel_without_checking(
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
    if (mt->remove_voxel_if_exists(coord))
    {
        (*itr) = mt->location();
        mt->location()->add_voxel_without_checking(
            particle_info_type(coord, ParticleID()));
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

bool LatticeSpaceVectorImpl::move_private(const private_coordinate_type& src,
        const private_coordinate_type& dest, const std::size_t candidate)
{
    return move_(src, dest, candidate).second;
}

bool LatticeSpaceVectorImpl::can_move(const private_coordinate_type& src,
        const private_coordinate_type& dest) const
{
    if (src == dest)
        return false;

    const MolecularTypeBase* src_mt(voxels_.at(src));
    if (src_mt->is_vacant())
        return false;

    MolecularTypeBase* dest_mt(voxels_.at(dest));

    if (dest_mt == border_)
        return false;

    if (dest_mt == periodic_)
        dest_mt = voxels_.at(apply_boundary_(dest));

    return (dest_mt == src_mt->location());
}

std::pair<LatticeSpaceVectorImpl::coordinate_type, bool>
    LatticeSpaceVectorImpl::move_to_neighbor(
        private_coordinate_type coord, Integer nrand)
{
    const private_coordinate_type neighbor(get_neighbor_private(coord, nrand));
    return move_(coord, neighbor);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool>
    LatticeSpaceVectorImpl::move_to_neighbor(
        particle_info_type& info, Integer nrand)
{
    const coordinate_type neighbor(get_neighbor_private(info.first, nrand));
    return move_(info, neighbor);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool>
    LatticeSpaceVectorImpl::move_(
        private_coordinate_type private_from, private_coordinate_type private_to,
        const std::size_t candidate)
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

    from_mt->replace_voxel(private_from, private_to, candidate);
    voxel_container::iterator from_itr(voxels_.begin() + private_from);
    (*from_itr) = to_mt;

    // to_mt->replace_voxel(private_to, particle_info_type(private_from, ParticleID()));
    to_mt->replace_voxel(private_to, private_from);
    voxel_container::iterator to_itr(voxels_.begin() + private_to);
    (*to_itr) = from_mt;

    return std::pair<private_coordinate_type, bool>(private_to, true);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool>
    LatticeSpaceVectorImpl::move_(
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

    info.first = private_to;
    voxel_container::iterator from_itr(voxels_.begin() + private_from);
    (*from_itr) = to_mt;

    // to_mt->replace_voxel(private_to, particle_info_type(private_from, ParticleID()));
    to_mt->replace_voxel(private_to, private_from);
    voxel_container::iterator to_itr(voxels_.begin() + private_to);
    (*to_itr) = from_mt;

    return std::pair<private_coordinate_type, bool>(private_to, true);
}

std::pair<LatticeSpaceVectorImpl::private_coordinate_type, bool>
    LatticeSpaceVectorImpl::move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand)
{
    const private_coordinate_type private_from(info.first);
    private_coordinate_type private_to(get_neighbor_private(private_from, nrand));

    //XXX: assert(private_from != private_to);
    //XXX: assert(from_mt == voxels_[private_from]);
    //XXX: assert(from_mt != vacant_);

    MolecularTypeBase* to_mt(voxels_[private_to]);

    if (to_mt != loc)
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

    to_mt->replace_voxel(private_to, private_from);
    // if (to_mt != vacant_) // (!to_mt->is_vacant())
    // {
    //     to_mt->replace_voxel(
    //         private_to, particle_info_type(private_from, ParticleID()));
    // }
    return std::make_pair(private_to, true);
}

const Particle LatticeSpaceVectorImpl::particle_at_private(
    private_coordinate_type coord) const
{
    const MolecularTypeBase* mt(voxels_.at(coord));
    return Particle(
        mt->species(),
        coordinate2position(private2coord(coord)),
        mt->radius(), mt->D());
}

Integer LatticeSpaceVectorImpl::num_voxels_exact(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        return 0;
    }
    const boost::shared_ptr<MolecularType>& mt((*itr).second);
    if (!mt->with_voxels())
    {
        // throw NotSupported(
        //     "num_voxels_exact for MolecularType with no voxel"
        //     " is not supporeted now");
        return count_voxels(mt);
    }
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
            const boost::shared_ptr<MolecularType>& mt((*itr).second);
            if (!mt->with_voxels())
            {
                // throw NotSupported(
                //     "num_voxels_exact for MolecularType with no voxel"
                //     " is not supporeted now");
                count += count_voxels(mt);
            }
            else
            {
                count += mt->size();
            }
        }
    }
    return count;
}

Integer LatticeSpaceVectorImpl::num_voxels() const
{
    Integer retval(0);
    for (spmap::const_iterator itr(spmap_.begin()); itr != spmap_.end(); ++itr)
    {
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        if (mt->with_voxels())
        {
            retval += (*itr).second->size();
        }
        else
        {
            retval += count_voxels(mt);  //XXX: too slow
        }
    }
    return retval;
}

/**
 * Change the Species at v.coordinate() to v.species.
 * The ParticleID must be kept after this update.
 */
void LatticeSpaceVectorImpl::update_voxel_private(const Voxel& v)
{
    const private_coordinate_type coord(v.coordinate());
    MolecularTypeBase* src_mt(voxels_.at(coord));
    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo

    if (src_mt->with_voxels() != new_mt->with_voxels())
    {
        throw NotSupported("ParticleID is needed/lost.");
    }

    new_mt->add_voxel_without_checking(src_mt->pop(coord));
    voxel_container::iterator itr(voxels_.begin() + coord);
    (*itr) = new_mt;

    // const LatticeSpaceVectorImpl::private_coordinate_type coord(v.coordinate());
    // MolecularTypeBase* src_mt(voxels_.at(coord));
    // // if (src_mt->is_vacant())
    // if (!src_mt->with_voxels())
    // {
    //     return false; // has no ParticleID. Call update_voxel_private with a ParticleID.
    // }

    // MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    // // if (new_mt->is_vacant())
    // if (!new_mt->with_voxels())
    // {
    //     return false; // has no ParticleID. Call remove_voxel.
    // }

    // // MolecularTypeBase::particle_info_type info(*(src_mt->find(coord)));
    // // src_mt->remove_voxel_if_exists(coord);
    // // new_mt->addVoxel(info);
    // // MolecularTypeBase::iterator position(src_mt->find(coord));
    // // new_mt->add_voxel_without_checking(*position);
    // // src_mt->remove_voxel(position);
    // new_mt->add_voxel_without_checking(src_mt->pop(coord));
    // voxel_container::iterator itr(voxels_.begin() + coord);
    // (*itr) = new_mt;
    // return true;
}

/*
 * Change the Species and coordinate of a Voxel with ParticleID, pid, to
 * v.species() and v.coordinate() respectively and return false.
 * If no Voxel with pid is found, create a new Voxel at v.coordiante() and return ture.
 */
bool LatticeSpaceVectorImpl::update_voxel_private(const ParticleID& pid, const Voxel& v)
{
    const LatticeSpaceVectorImpl::private_coordinate_type& to_coord(v.coordinate());
    if (!is_in_range_private(to_coord))
    {
        throw NotSupported("Out of bounds");
    }

    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    MolecularTypeBase* dest_mt(get_molecular_type(to_coord));

    if (dest_mt != new_mt->location())
    {
        throw NotSupported(
            "Mismatch in the location. Failed to place '"
            + new_mt->species().serial() + "' to '"
            + dest_mt->species().serial() + "'.");
    }

    const LatticeSpaceVectorImpl::private_coordinate_type
        from_coord(pid != ParticleID() ? get_coord(pid) : -1);
    if (from_coord != -1)
    {
        // move
        MolecularTypeBase* src_mt(voxels_.at(from_coord));
        src_mt->remove_voxel_if_exists(from_coord);

        //XXX: use location?
        dest_mt->replace_voxel(to_coord, from_coord);
        voxel_container::iterator from_itr(voxels_.begin() + from_coord);
        (*from_itr) = dest_mt;

        new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
        voxel_container::iterator to_itr(voxels_.begin() + to_coord);
        (*to_itr) = new_mt;
        return false;
    }

    // new
    dest_mt->remove_voxel_if_exists(to_coord);

    new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
    voxel_container::iterator to_itr(voxels_.begin() + to_coord);
    (*to_itr) = new_mt;
    return true;

    // const LatticeSpaceVectorImpl::private_coordinate_type& to_coord(v.coordinate());
    // if (!is_in_range_private(to_coord))
    // {
    //     return false;
    // }

    // MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    // if (new_mt->is_vacant())
    // {
    //     return false; // has no ParticleID.
    // }

    // MolecularTypeBase* dest_mt(get_molecular_type(to_coord));
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
    // }

    // if (pid != ParticleID())
    // {
    //     if (!new_mt->with_voxels())
    //     {
    //         return false; // has no ParticleID
    //     }

    //     const LatticeSpaceVectorImpl::private_coordinate_type from_coord(get_coord(pid));
    //     if (from_coord != -1)
    //     {
    //         MolecularTypeBase* src_mt(voxels_.at(from_coord));
    //         // assert(src_mt->with_voxels());
    //         src_mt->remove_voxel_if_exists(from_coord);
    //         voxel_container::iterator itr(voxels_.begin() + from_coord);
    //         (*itr) = dest_mt;
    //         dest_mt->replace_voxel(to_coord, from_coord);
    //         new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
    //         itr = voxels_.begin() + to_coord;
    //         (*itr) = new_mt;
    //         return true;
    //     }
    // }

    // new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
    // voxel_container::iterator to_itr(voxels_.begin() + to_coord);
    // (*to_itr) = new_mt;
    // dest_mt->remove_voxel_if_exists(to_coord);
    // return true;
}

bool LatticeSpaceVectorImpl::update_voxel_private_without_checking(const ParticleID& pid, const Voxel& v)
{
    const LatticeSpaceVectorImpl::private_coordinate_type& to_coord(v.coordinate());
    if (!is_in_range_private(to_coord))
    {
        throw NotSupported("Out of bounds");
    }

    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    MolecularTypeBase* dest_mt(get_molecular_type(to_coord));

    const LatticeSpaceVectorImpl::private_coordinate_type
        from_coord(pid != ParticleID() ? get_coord(pid) : -1);
    if (from_coord != -1)
    {
        // move
        MolecularTypeBase* src_mt(voxels_.at(from_coord));
        src_mt->remove_voxel_if_exists(from_coord);

        //XXX: use location?
        dest_mt->replace_voxel(to_coord, from_coord);
        voxel_container::iterator from_itr(voxels_.begin() + from_coord);
        (*from_itr) = dest_mt;

        new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
        voxel_container::iterator to_itr(voxels_.begin() + to_coord);
        (*to_itr) = new_mt;
        return false;
    }

    // new
    dest_mt->remove_voxel_if_exists(to_coord);

    new_mt->add_voxel_without_checking(particle_info_type(to_coord, pid));
    voxel_container::iterator to_itr(voxels_.begin() + to_coord);
    (*to_itr) = new_mt;
    return true;

}

bool LatticeSpaceVectorImpl::make_structure_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    spmap::iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
    {
        return false;
    }

    MolecularTypeBase* location;
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
            // XXX: A MolecularTypeBase for the structure (location) must be allocated
            // XXX: before the allocation of a Species on the structure.
            // XXX: The MolecularTypeBase cannot be automatically allocated at the time
            // XXX: because its MoleculeInfo is unknown.
            // XXX: LatticeSpaceVectorImpl::load will raise a problem about this issue.
            // XXX: In this implementation, the MolecularTypeBase for a structure is
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
    std::pair<spmap::iterator, bool> retval(spmap_.insert(std::make_pair(sp, mt)));
    return retval.second;
}

bool LatticeSpaceVectorImpl::make_molecular_type(const Species& sp, Real radius, Real D, const std::string loc)
{
    spmap::iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
        return false;
    std::pair<spmap::iterator, bool> retval(__get_molecular_type(Voxel(sp, 0, radius, D, loc)));
    return retval.second;
}

bool LatticeSpaceVectorImpl::add_voxels(const Species sp, std::vector<std::pair<ParticleID, coordinate_type> > voxels)
{
    // this function doesn't check location.
    MolecularTypeBase *mtb;
    try
    {
        mtb = find_molecular_type(sp);
    }
    catch (NotFound &e)
    {
        return false;
    }
    for (std::vector<std::pair<ParticleID, coordinate_type> >::iterator itr(voxels.begin());
            itr != voxels.end(); ++itr)
    {
        const ParticleID pid((*itr).first);
        const private_coordinate_type coord(coord2private((*itr).second));
        MolecularTypeBase* src_mt(get_molecular_type(coord));
        src_mt->remove_voxel_if_exists(coord);
        mtb->add_voxel_without_checking(particle_info_type(coord, pid));
        voxel_container::iterator vitr(voxels_.begin() + coord);
        (*vitr) = mtb;
    }
    return true;
}

Integer LatticeSpaceVectorImpl::count_voxels(
    const boost::shared_ptr<MolecularType>& mt) const
{
    return static_cast<Integer>(
        std::count(voxels_.begin(), voxels_.end(), mt.get()));
}

} // ecell4
