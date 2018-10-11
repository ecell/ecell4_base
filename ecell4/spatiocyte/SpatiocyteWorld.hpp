#ifndef ECELL4_LATTICE_LATTICE_WORLD_HPP
#define ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <sstream>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/LatticeSpaceCellListImpl.hpp>
#include <ecell4/core/LatticeSpaceVectorImpl.hpp>
#include <ecell4/core/OffLatticeSpace.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/WorldInterface.hpp>

#include "OneToManyMap.hpp"
#include "Voxel.hpp"
#include "SpatiocyteReactions.hpp"

namespace ecell4
{

namespace spatiocyte
{

struct MoleculeInfo
{
    const Real radius;
    const Real D;
    const std::string loc;
};

class SpatiocyteWorld
    : public WorldInterface
{
public:

    typedef LatticeSpaceVectorImpl default_root_type;

    typedef MoleculeInfo molecule_info_type;

    typedef VoxelSpaceBase::coordinate_id_pair_type coordinate_id_pair_type;
    typedef VoxelSpaceBase::coordinate_type coordinate_type;

    typedef boost::shared_ptr<VoxelSpaceBase> space_type;
    typedef std::vector<space_type> space_container_type;

public:

    /*
     * Constructors
     */
    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
        : rng_(rng)
    {
        spaces_.push_back(space_type(new default_root_type(edge_lengths, voxel_radius)));
        size_ = get_root()->size();
    }

    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius)
    {
        spaces_.push_back(space_type(new default_root_type(edge_lengths, voxel_radius)));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
        size_ = get_root()->size();
    }

    SpatiocyteWorld(const Real3& edge_lengths = Real3(1, 1, 1))
    {
        // XXX: sloppy default
        spaces_.push_back(space_type(new default_root_type(edge_lengths, edge_lengths[0] / 100)));
        size_ = get_root()->size();
        rng_ = boost::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const std::string filename)
    {
        // XXX: sloppy default
        spaces_.push_back(space_type(new default_root_type(Real3(1, 1, 1), 1 / 100)));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }

    SpatiocyteWorld(VoxelSpaceBase* space, const boost::shared_ptr<RandomNumberGenerator>& rng)
        : rng_(rng)
    {
        spaces_.push_back(space_type(space));
        size_ = get_root()->size();
    }

    void add_space(VoxelSpaceBase *space);

    const Real t() const
    {
        Real time(0.0);

        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            time = std::max(time, (*itr)->t());
        }

        return time;
    }

    void set_t(const Real& t)
    {
        for (space_container_type::iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            (*itr)->set_t(t);
        }
    }

    void save(const std::string& filename) const
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File> fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        sidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group> group(new H5::Group(fout->createGroup("LatticeSpace")));
        get_root()->save_hdf5(group.get()); // TODO
        extras::save_version_information(fout.get(), std::string("ecell4-spatiocyte-") + std::string(ECELL4_VERSION));
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
    }

    void load(const std::string& filename)
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));

        const std::string required = "ecell4-spatiocyte-4.1.0";
        try
        {
            const std::string version = extras::load_version_information(*fin);
            if (!extras::check_version_information(version, required))
            {
                std::stringstream ss;
                ss << "The version of the given file [" << version
                    << "] is too old. [" << required << "] or later is required.";
                throw NotSupported(ss.str());
            }
        }
        catch(H5::GroupIException not_found_error)
        {
            throw NotFound("No version information was found.");
        }

        const H5::Group group(fin->openGroup("LatticeSpace"));
        get_root()->load_hdf5(group); // TODO
        sidgen_.load(*fin);
        rng_->load(*fin);
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
    }

    // Integer num_species() const
    // {
    //     Integer total(0);
    //     for (space_container_type::const_iterator itr(spaces_.begin());
    //          itr != spaces_.end(); ++itr)
    //     {
    //         total += (*itr)->num_species();
    //     }
    //     return total;
    // }

    bool has_species(const Species &sp) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_species(sp))
                return true;
        }
        return false;
    }

    Integer num_molecules(const Species& sp) const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_molecules(sp);
        }
        return total;
    }

    Integer num_molecules_exact(const Species& sp) const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_molecules_exact(sp);
        }
        return total;
    }

    Real get_value(const Species& sp) const
    {
        Real value(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            value += (*itr)->get_value(sp);
        }
        return value;
    }

    Real get_value_exact(const Species& sp) const
    {
        Real value(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            value += (*itr)->get_value_exact(sp);
        }
        return value;
    }

    const Real3& edge_lengths() const
    {
        return get_root()->edge_lengths();
    }

    Integer num_particles() const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_particles();
        }
        return total;
    }

    Integer num_particles(const Species& sp) const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_particles(sp);
        }
        return total;
    }

    Integer num_particles_exact(const Species& sp) const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_particles_exact(sp);
        }
        return total;
    }

    bool has_particle(const ParticleID& pid) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_particle(pid))
                return true;
        }
        return false;
    }

    // Suggests: Rename to 'find_particle'
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_particle(pid))
                return (*itr)->get_particle(pid);
        }
        throw "No particle corresponding to a given ParticleID is found.";
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles() const
    {
        std::vector<std::pair<ParticleID, Particle> > list;
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<std::pair<ParticleID, Particle> > particles((*itr)->list_particles());
            list.insert(list.end(), particles.begin(), particles.end());
        }
        return list;
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, Particle> > list;
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<std::pair<ParticleID, Particle> > particles((*itr)->list_particles(sp));
            list.insert(list.end(), particles.begin(), particles.end());
        }

        return list;
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles_exact(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, Particle> > list;
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<std::pair<ParticleID, Particle> > particles(
                    (*itr)->list_particles_exact(sp));
            list.insert(list.end(), particles.begin(), particles.end());
        }
        return list;
    }

    std::vector<std::pair<ParticleID, Particle> > list_structure_particles() const;
    std::vector<std::pair<ParticleID, Particle> > list_non_structure_particles() const;

    Real voxel_radius() const
    {
        return get_root()->voxel_radius();
    }

    Real voxel_volume() const
    {
        return get_root()->voxel_volume();
    }

    Real get_volume(const Species& sp) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_species(sp) && (*itr)->find_molecule_pool(sp)->is_structure())
                return (*itr)->get_volume(sp);
        }
        return 0.0;
    }

    const Real volume() const
    {
        return get_root()->volume();
    }

    Real unit_area() const
    {
        return get_root()->unit_area();
    }

    // TODO
    boost::shared_ptr<VoxelPool> vacant() const {
        return get_root()->vacant();
    }

    bool has_voxel(const ParticleID& pid) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_voxel(pid))
                return true;
        }
        return false;
    }

    Integer num_voxels() const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_voxels();
        }
        return total;
    }

    Integer num_voxels(const Species& sp) const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_voxels(sp);
        }
        return total;
    }

    Integer num_voxels_exact(const Species& sp) const
    {
        Integer total(0);
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            total += (*itr)->num_voxels_exact(sp);
        }
        return total;
    }

    std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels() const
    {
        std::vector<std::pair<ParticleID, ParticleVoxel> > list;
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<std::pair<ParticleID, ParticleVoxel> > voxels((*itr)->list_voxels());
            list.insert(list.end(), voxels.begin(), voxels.end());
        }
        return list;
    }

    std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, ParticleVoxel> > list;
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<std::pair<ParticleID, ParticleVoxel> > voxels((*itr)->list_voxels(sp));
            list.insert(list.end(), voxels.begin(), voxels.end());
        }
        return list;
    }

    std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels_exact(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, ParticleVoxel> > list;
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<std::pair<ParticleID, ParticleVoxel> > voxels(
                    (*itr)->list_voxels_exact(sp));
            list.insert(list.end(), voxels.begin(), voxels.end());
        }
        return list;
    }

    boost::optional<ParticleVoxel> find_voxel(const ParticleID& pid) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if (boost::optional<ParticleVoxel> voxel = (*itr)->find_voxel(pid))
            {
                return voxel;
            }
        }

        return boost::none;
    }

    std::pair<ParticleID, Species> get_voxel_at(const Voxel& voxel) const
    {
        std::pair<ParticleID, ParticleVoxel> id_voxel_pair(
                voxel.space.lock()->get_voxel_at(voxel.coordinate));
        return std::make_pair(id_voxel_pair.first, id_voxel_pair.second.species);
    }

    boost::shared_ptr<VoxelPool> find_voxel_pool(const Species& species)
    {
        for (space_container_type::iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_species(species))
                return (*itr)->find_voxel_pool(species);
        }
        // create VoxelPool TODO
        return get_root()->find_voxel_pool(species);
        // throw "No VoxelPool corresponding to a given Species is found";
    }

    boost::shared_ptr<const VoxelPool> find_voxel_pool(const Species& species) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_species(species))
                return (*itr)->find_voxel_pool(species);
        }
        throw "No VoxelPool corresponding to a given Species is found";
    }

    bool has_molecule_pool(const Species& species) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_molecule_pool(species))
                return true;
        }
        return false;
    }

    boost::shared_ptr<MoleculePool> find_molecule_pool(const Species& species)
    {
        for (space_container_type::iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_molecule_pool(species))
                return (*itr)->find_molecule_pool(species);
        }
        throw "No MoleculePool corresponding to a given Species is found";
    }

    boost::shared_ptr<const MoleculePool> find_molecule_pool(const Species& species) const
    {
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_molecule_pool(species))
                return (*itr)->find_molecule_pool(species);
        }
        throw "No MoleculePool corresponding to a given Species is found";
    }

    /*
     * Coordinate Transformation
     */
    Voxel position2voxel(const Real3& pos) const
    {
        return Voxel(get_root(), get_root()->position2coordinate(pos));
    }

    /*
     * ParticleVoxel Manipulation
     */
    bool update_voxel(const ParticleID& pid, ParticleVoxel v)
    {
        for (space_container_type::iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_voxel(pid))
                return (*itr)->update_voxel(pid, v);
        }

        return get_space(v.coordinate)->update_voxel(pid, v);
    }

    bool remove_voxel(const ParticleID& pid)
    {
        for (space_container_type::iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_voxel(pid))
                return (*itr)->remove_voxel(pid);
        }
        return false;
    }

    // Deprecated
    bool can_move(const Voxel& src, const Voxel& dst) const
    {
        if (!(src.space < dst.space) && !(dst.space < src.space))
            return src.space.lock()->can_move(src.coordinate, dst.coordinate);

        return false;
    }

    // Deprecated
    bool move(const Voxel& src, const Voxel& dst, const std::size_t candidate=0)
    {
        if (!(src.space < dst.space) && !(dst.space < src.space))
            return src.space.lock()->move(src.coordinate, dst.coordinate, candidate);

        return false;

    }

    const Integer size() const
    {
        return size_;
    }

    const Integer3 shape() const
    {
        return get_root()->shape();
    }

    /*
     * SpatiocyteWorld API
     */

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    MoleculeInfo get_molecule_info(const Species& sp) const
    {
        const bool with_D(sp.has_attribute("D"));
        const bool with_radius(sp.has_attribute("radius"));
        const bool with_loc(sp.has_attribute("location"));

        Real radius(voxel_radius()), D(0.0);
        std::string loc("");

        if (with_D)
        {
            D = sp.get_attribute_as<Real>("D");
        }

        if (with_radius)
        {
            radius = sp.get_attribute_as<Real>("radius");
        }

        if (with_loc)
        {
            loc = sp.get_attribute_as<std::string>("location");
        }

        if (!(with_D && with_radius))  //XXX: with_loc?
        {
            if (boost::shared_ptr<Model> bound_model = lock_model())
            {
                Species newsp(bound_model->apply_species_attributes(sp));

                if (!with_D && newsp.has_attribute("D"))
                {
                    D = newsp.get_attribute_as<Real>("D");
                }

                if (!with_radius && newsp.has_attribute("radius"))
                {
                    radius = newsp.get_attribute_as<Real>("radius");
                }

                if (!with_loc && newsp.has_attribute("location"))
                {
                    loc = newsp.get_attribute_as<std::string>("location");
                }
            }
        }

        MoleculeInfo info = {radius, D, loc};
        return info;
    }

    // bool has_species_exact(const Species &sp) const
    // {
    //     return get_root()->has_species_exact(sp);
    // }

    void set_value(const Species& sp, const Real value);


    /**
     * create and add a new particle
     * @param p a particle
     * @return a pair of a pair of pid (a particle id) and p (a particle)
     * and bool (if it's succeeded or not)
     */
    boost::optional<ParticleID>
    new_particle(const Particle& p)
    {
        // ParticleID pid(sidgen_());
        // const bool is_succeeded(update_particle(pid, p));
        // return std::make_pair(get_particle(pid), is_succeeded);
        const molecule_info_type minfo(get_molecule_info(p.species()));
        const Voxel voxel(position2voxel(p.position()));

        if (voxel.get_voxel_pool()->species().serial() != minfo.loc)
            return boost::none;

        if (boost::optional<ParticleID> pid = new_voxel(p.species(), voxel))
            return *pid;

        return boost::none;
    }

    boost::optional<ParticleID>
    new_particle(const Species& sp, const Real3& pos)
    {
        const MoleculeInfo info(get_molecule_info(sp));
        return new_particle(Particle(sp, pos, info.radius, info.D));
    }

    bool remove_particle(const ParticleID& pid)
    {
        for (space_container_type::iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr)->has_particle(pid))
                return (*itr)->remove_particle(pid);
        }
        return false;
    }

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        const molecule_info_type minfo(get_molecule_info(p.species()));
        return update_voxel(pid,
                ParticleVoxel(p.species(),
                              position2voxel(p.position()).coordinate,
                              p.radius(), p.D(), minfo.loc));
    }

    std::vector<Species> list_species() const
    {
        std::vector<Species> list;
        for (space_container_type::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<Species> species((*itr)->list_species());
            list.insert(list.end(), species.begin(), species.end());
        }
        return list;
    }

    std::vector<Species> list_non_structure_species() const;
    std::vector<Species> list_structure_species() const;
    // std::vector<coordinate_type> list_coords(const Species& sp) const;

    boost::optional<ParticleID> new_voxel(const Species& sp, const Voxel& voxel)
    {
        boost::shared_ptr<VoxelSpaceBase> space(voxel.space.lock());
        if (!space->has_species(sp))
        {
            const molecule_info_type minfo(get_molecule_info(sp));
            space->make_molecular_type(sp, minfo.radius, minfo.D, minfo.loc);
        }

        ParticleID pid(sidgen_());

        if (space->add_voxel(sp, pid, voxel.coordinate))
            return pid;

        return boost::none;
    }

    boost::optional<ParticleID>
    new_voxel_structure(const Species& sp, const Voxel& voxel)
    {
        boost::shared_ptr<VoxelSpaceBase> space(voxel.space.lock());
        if (!space->has_species(sp))
        {
            const molecule_info_type minfo(get_molecule_info(sp));
            space->make_molecular_type(sp, minfo.radius, minfo.D, minfo.loc);
        }

        ParticleID pid;

        if (space->add_voxel(sp, pid, voxel.coordinate))
            return pid;

        return boost::none;
    }

    bool add_molecules(const Species& sp, const Integer& num);
    bool add_molecules(const Species& sp, const Integer& num, const boost::shared_ptr<const Shape> shape);
    Integer add_structure(const Species& sp, const boost::shared_ptr<const Shape> shape);

    void remove_molecules(const Species& sp, const Integer& num);
    // void remove_molecules_exact(const Species& sp, const Integer& num);

    boost::optional<Voxel>
    check_neighbor(const Voxel& voxel, const std::string& loc);

    // bool update_molecule(coordinate_type at, Species species);

    const Species& draw_species(const Species& pttrn) const;

    boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    const molecule_info_type get_molecule_info(boost::shared_ptr<const VoxelPool> mt) const
    {
        const molecule_info_type info = {mt->radius(), mt->D(), get_location_serial(mt)};
        return info;
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to SpatiocyteWorld"
                    << std::endl;
            }
        }

        model_ = model;
    }

    boost::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
    }

    Shape::dimension_kind get_dimension(const Species& species) const
    {
        return extras::get_dimension_from_model(species, lock_model());
    }

    /**
     * static members
     */
    static inline Real calculate_voxel_volume(const Real r)
    {
        return VoxelSpaceBase::calculate_voxel_volume(r);
    }

    static inline Real3 calculate_hcp_lengths(const Real voxel_radius)
    {
        return VoxelSpaceBase::calculate_hcp_lengths(voxel_radius);
    }

    static inline Integer3 calculate_shape(const Real3& edge_lengths, const Real& voxel_radius)
    {
        return VoxelSpaceBase::calculate_shape(edge_lengths, voxel_radius, true);
    }

    static inline Real calculate_volume(const Real3& edge_lengths, const Real& voxel_radius)
    {
        return VoxelSpaceBase::calculate_volume(edge_lengths, voxel_radius, true);
    }

protected:

    bool is_inside(const coordinate_type& coord) const
    {
        return get_root()->is_inside(coord);
    }

    space_type get_root() const
    {
        return spaces_.at(0);
    }

    space_type get_space(coordinate_type& coordinate)
    {
        for (space_container_type::iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if (coordinate < (*itr)->size())
                return *itr;

            coordinate -= (*itr)->size();
        }
        return space_type();
    }

    Integer add_structure2(const Species& sp, const boost::shared_ptr<const Shape> shape);
    Integer add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape);
    bool is_surface_voxel(const Voxel& voxel, const boost::shared_ptr<const Shape> shape) const;

public:

    // TODO: Calling this function is invalid, and this should be removed.
    Voxel coordinate2voxel(const coordinate_type& coordinate)
    {
        coordinate_type coord(coordinate);
        return Voxel(get_space(coord), coord);
    }

protected:

    std::size_t size_;
    space_container_type spaces_;

    OneToManyMap<coordinate_type> interfaces_;
    OneToManyMap<coordinate_type> neighbors_;

    boost::shared_ptr<RandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> sidgen_;

    boost::weak_ptr<Model> model_;
};

inline
SpatiocyteWorld*
create_spatiocyte_world_cell_list_impl(
        const Real3& edge_lengths,
        const Real& voxel_radius,
        const Integer3& matrix_sizes,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new SpatiocyteWorld(
        new LatticeSpaceCellListImpl(edge_lengths, voxel_radius, matrix_sizes), rng);
}

inline
SpatiocyteWorld*
create_spatiocyte_world_vector_impl(
        const Real3& edge_lengths,
        const Real& voxel_radius,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new SpatiocyteWorld(
        new LatticeSpaceVectorImpl(edge_lengths, voxel_radius), rng);
}

inline
SpatiocyteWorld*
allocate_spatiocyte_world_square_offlattice_impl(
        const Real edge_length,
        const Real& voxel_radius,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    OffLatticeSpace::position_container positions;
    OffLatticeSpace::coordinate_pair_list_type adjoining_pairs;

    // const std::size_t num_row(int((edge_length - (2 + sqrt(3)) * voxel_radius) /
    //                          (2 * sqrt(3) * voxel_radius)) + 1);
    const std::size_t num_row(int(edge_length / (2 * sqrt(3) * voxel_radius)));
    const std::size_t num_col(int(edge_length / (2 * voxel_radius)));

    for (std::size_t row(0); row < num_row; ++row)
        for (std::size_t col(0); col < num_col; ++col)
        {
            // 2 * (row * num_col + col)
            positions.push_back(
                    Real3(2*col, 2*row*sqrt(3), 0) * voxel_radius);

            // 2 * (row * num_col + col + 1)
            positions.push_back(
                    Real3(2*col, (2*row+1)*sqrt(3), 0) * voxel_radius);

            const int index(2 * (row * num_col + col));

            const std::size_t next_col((col + 1) % num_col);
            const std::size_t next_row((row + 1) % num_row);

            const int right(2 * (row * num_col + next_col));
            const int bottom(2 * (next_row * num_col + col));
            const int right_bottom(2 * (next_row * num_col + next_col));

            adjoining_pairs.push_back(std::make_pair(index,   index+1));
            adjoining_pairs.push_back(std::make_pair(index,   right));
            adjoining_pairs.push_back(std::make_pair(index+1, right));
            adjoining_pairs.push_back(std::make_pair(index+1, right+1));
            adjoining_pairs.push_back(std::make_pair(index+1, bottom));
            adjoining_pairs.push_back(std::make_pair(index+1, right_bottom));
        }

    OffLatticeSpace *space = new OffLatticeSpace(voxel_radius, positions, adjoining_pairs, Shape::TWO);
    space->set_lengths(Real3(2*num_col, 2*sqrt(3)*num_row, 2) * voxel_radius);

    return new SpatiocyteWorld(space, rng);
}

/**
 * Alias functions for Cython
 */

inline SpatiocyteWorld* create_spatiocyte_world_cell_list_impl_alias(
    const Real3& edge_lengths, const Real& voxel_radius,
    const Integer3& matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return create_spatiocyte_world_cell_list_impl(
        edge_lengths, voxel_radius, matrix_sizes, rng);
}

inline SpatiocyteWorld* create_spatiocyte_world_vector_impl_alias(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return create_spatiocyte_world_vector_impl(edge_lengths, voxel_radius, rng);
}

} // spatiocyte

} // ecell4

#endif /* ECELL4_LATTICE_LATTICE_WORLD_HPP */
