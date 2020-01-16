#ifndef ECELL4_LATTICE_LATTICE_WORLD_HPP
#define ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include <ecell4/core/LatticeSpaceCellListImpl.hpp>
#include <ecell4/core/LatticeSpaceVectorImpl.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/OffLatticeSpace.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/WorldInterface.hpp>
#include <ecell4/core/extras.hpp>

#include "OneToManyMap.hpp"
#include "SpatiocyteReactions.hpp"
#include "Voxel.hpp"

namespace ecell4
{

namespace spatiocyte
{

struct MoleculeInfo
{
    const Real radius;
    const Real D;
    const std::string loc;
    const Shape::dimension_kind dimension;
};

class SpatiocyteWorld : public WorldInterface
{
public:
    typedef LatticeSpaceVectorImpl default_root_type;

    typedef VoxelSpaceBase::coordinate_id_pair_type coordinate_id_pair_type;
    typedef VoxelSpaceBase::coordinate_type coordinate_type;

    typedef boost::shared_ptr<VoxelSpaceBase> space_type;
    typedef std::vector<space_type> space_container_type;

public:
    /*
     * Constructors
     */
    SpatiocyteWorld(const Real3 &edge_lengths, const Real &voxel_radius,
                    const boost::shared_ptr<RandomNumberGenerator> &rng)
        : rng_(rng)
    {
        spaces_.push_back(
            space_type(new default_root_type(edge_lengths, voxel_radius)));
        size_ = get_root()->size();
    }

    SpatiocyteWorld(const Real3 &edge_lengths, const Real &voxel_radius)
    {
        spaces_.push_back(
            space_type(new default_root_type(edge_lengths, voxel_radius)));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
        size_ = get_root()->size();
    }

    SpatiocyteWorld(const Real3 &edge_lengths = Real3(1, 1, 1))
    {
        // XXX: sloppy default
        spaces_.push_back(space_type(
            new default_root_type(edge_lengths, edge_lengths[0] / 100)));
        size_ = get_root()->size();
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const std::string filename)
    {
        // XXX: sloppy default
        spaces_.push_back(
            space_type(new default_root_type(Real3(1, 1, 1), 1 / 100)));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }

    SpatiocyteWorld(VoxelSpaceBase *space,
                    const boost::shared_ptr<RandomNumberGenerator> &rng)
        : rng_(rng)
    {
        spaces_.push_back(space_type(space));
        size_ = get_root()->size();
    }

    void add_space(std::unique_ptr<VoxelSpaceBase> space);

    const Real t() const
    {
        Real time(0.0);

        for (const auto &space : spaces_)
        {
            time = std::max(time, space->t());
        }

        return time;
    }

    void set_t(const Real &t)
    {
        for (auto &space : spaces_)
        {
            space->set_t(t);
        }
    }

    void save(const std::string &filename) const
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File> fout(
            new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        sidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group> group(
            new H5::Group(fout->createGroup("LatticeSpace")));
        get_root()->save_hdf5(group.get()); // TODO
        extras::save_version_information(fout.get(),
                                         std::string("ecell4-spatiocyte-") +
                                             std::string(VERSION_INFO));
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
    }

    void load(const std::string &filename)
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File> fin(
            new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));

        const std::string required = "ecell4-spatiocyte-0.0";
        try
        {
            const std::string version = extras::load_version_information(*fin);
            if (!extras::check_version_information(version, required))
            {
                std::stringstream ss;
                ss << "The version of the given file [" << version
                   << "] is too old. [" << required
                   << "] or later is required.";
                throw NotSupported(ss.str());
            }
        }
        catch (H5::GroupIException not_found_error)
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
        for (const auto &space : spaces_)
        {
            if (space->has_species(sp))
                return true;
        }
        return false;
    }

    Integer num_molecules(const Species &sp) const
    {
        Integer total(0);
        for (const auto &space : spaces_)
        {
            total += space->num_molecules(sp);
        }
        return total;
    }

    Integer num_molecules_exact(const Species &sp) const
    {
        Integer total(0);
        for (const auto &space : spaces_)
        {
            total += space->num_molecules_exact(sp);
        }
        return total;
    }

    Real get_value(const Species &sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    Real get_value_exact(const Species &sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

    const Real3 &edge_lengths() const { return get_root()->edge_lengths(); }

    Integer num_particles() const { return num_voxels(); }

    Integer num_particles(const Species &sp) const { return num_voxels(sp); }

    Integer num_particles_exact(const Species &sp) const
    {
        return num_voxels_exact(sp);
    }

    boost::optional<Voxel> get_voxel(const ParticleID &pid) const
    {
        for (const auto &space : spaces_)
        {
            if (const auto coordinate = space->get_coordinate(pid))
                return Voxel(space, *coordinate);
        }
        return boost::none;
    }

    bool has_particle(const ParticleID &pid) const { return has_voxel(pid); }

    // Suggests: Rename to 'find_particle'
    std::pair<ParticleID, Particle> get_particle(const ParticleID &pid) const
    {
        for (const auto &space : spaces_)
        {
            if (const auto &view = space->find_voxel(pid))
            {
                return std::make_pair(pid, gen_particle_from(space, *view));
            }
        }
        throw "No particle corresponding to a given ParticleID is found.";
    }

    std::vector<std::pair<ParticleID, Particle>> list_particles() const
    {
        return list_particles_private(
            [](const space_type &space) { return space->list_voxels(); });
    }

    std::vector<std::pair<ParticleID, Particle>>
    list_particles(const Species &sp) const
    {
        return list_particles_private(
            [&sp](const space_type &space) { return space->list_voxels(sp); });
    }

    std::vector<std::pair<ParticleID, Particle>>
    list_particles_exact(const Species &sp) const
    {
        return list_particles_private([&sp](const space_type &space) {
            return space->list_voxels_exact(sp);
        });
    }

    std::vector<std::pair<ParticleID, Particle>>
    list_structure_particles() const;
    std::vector<std::pair<ParticleID, Particle>>
    list_non_structure_particles() const;

    Real voxel_radius() const { return get_root()->voxel_radius(); }

    Real voxel_volume() const { return get_root()->voxel_volume(); }

    Real get_volume(const Species &sp) const
    {
        for (const auto &space : spaces_)
        {
            if (space->has_species(sp) &&
                space->find_molecule_pool(sp)->is_structure())
                return space->get_volume(sp);
        }
        return 0.0;
    }

    const Real volume() const { return get_root()->volume(); }

    Real unit_area() const { return get_root()->unit_area(); }

    // TODO
    boost::shared_ptr<VoxelPool> vacant() const { return get_root()->vacant(); }

    bool has_voxel(const ParticleID &pid) const
    {
        for (const auto &space : spaces_)
        {
            if (space->has_voxel(pid))
                return true;
        }
        return false;
    }

    Integer num_voxels() const
    {
        Integer total(0);
        for (const auto &space : spaces_)
        {
            total += space->num_voxels();
        }
        return total;
    }

    Integer num_voxels(const Species &sp) const
    {
        Integer total(0);
        for (const auto &space : spaces_)
        {
            total += space->num_voxels(sp);
        }
        return total;
    }

    Integer num_voxels_exact(const Species &sp) const
    {
        Integer total(0);
        for (const auto &space : spaces_)
        {
            total += space->num_voxels_exact(sp);
        }
        return total;
    }

    std::vector<ParticleBase<Voxel>> list_voxels() const
    {
        return list_voxels_private(
            [](const space_type &space) { return space->list_voxels(); });
    }

    std::vector<ParticleBase<Voxel>> list_voxels(const Species &sp) const
    {
        return list_voxels_private(
            [&sp](const space_type &space) { return space->list_voxels(sp); });
    }

    std::vector<ParticleBase<Voxel>> list_voxels_exact(const Species &sp) const
    {
        return list_voxels_private([&sp](const space_type &space) {
            return space->list_voxels_exact(sp);
        });
    }

    Species get_species_at(const Voxel &voxel) const
    {
        return voxel.get_voxel_pool()->species();
    }

    bool has_particle_at(const Voxel &voxel) const
    {
        return !voxel.get_voxel_pool()->is_vacant();
    }

    std::pair<ParticleID, Species> get_voxel_at(const Voxel &voxel) const
    {
        const auto view(voxel.space.lock()->get_voxel_at(voxel.coordinate));
        return std::make_pair(view.pid, view.species);
    }

    boost::shared_ptr<VoxelPool> find_voxel_pool(const Species &species)
    {
        for (const auto &space : spaces_)
        {
            if (space->has_species(species))
                return space->find_voxel_pool(species);
        }
        // create VoxelPool TODO
        return get_root()->find_voxel_pool(species);
        // throw "No VoxelPool corresponding to a given Species is found";
    }

    boost::shared_ptr<const VoxelPool>
    find_voxel_pool(const Species &species) const
    {
        for (const auto &space : spaces_)
        {
            if (space->has_species(species))
                return space->find_voxel_pool(species);
        }
        throw "No VoxelPool corresponding to a given Species is found";
    }

    boost::optional<std::pair<space_type, boost::shared_ptr<VoxelPool>>>
    find_space_and_voxel_pool(const Species &species) const
    {
        for (const auto &space : spaces_)
        {
            if (space->has_species(species))
            {
                const auto voxel_pool = space->find_voxel_pool(species);
                return std::pair<space_type, boost::shared_ptr<VoxelPool>>(
                    space, voxel_pool);
            }
        }
        return boost::none;
    }

    bool has_molecule_pool(const Species &species) const
    {
        for (const auto &space : spaces_)
        {
            if (space->has_molecule_pool(species))
                return true;
        }
        return false;
    }

    boost::shared_ptr<MoleculePool> find_molecule_pool(const Species &species)
    {
        for (const auto &space : spaces_)
        {
            if (space->has_molecule_pool(species))
                return space->find_molecule_pool(species);
        }
        throw "No MoleculePool corresponding to a given Species is found";
    }

    boost::shared_ptr<const MoleculePool>
    find_molecule_pool(const Species &species) const
    {
        for (const auto &space : spaces_)
        {
            if (space->has_molecule_pool(species))
                return space->find_molecule_pool(species);
        }
        throw "No MoleculePool corresponding to a given Species is found";
    }

    boost::optional<std::pair<space_type, boost::shared_ptr<MoleculePool>>>
    find_space_and_molecule_pool(const Species &species) const
    {
        for (const auto &space : spaces_)
        {
            if (space->has_molecule_pool(species))
            {
                const auto molecule_pool(space->find_molecule_pool(species));
                return std::pair<space_type, boost::shared_ptr<MoleculePool>>(
                    space, molecule_pool);
            }
        }
        return boost::none;
    }

    /*
     * Coordinate Transformation
     */
    Voxel get_voxel_nearby(const Real3 &pos) const
    {
        return Voxel(get_root(), get_root()->position2coordinate(pos));
    }

    /*
     * Voxel Manipulation
     */
    bool update_voxel(const ParticleID &pid, const Species species,
                      const Voxel voxel)
    {
        const MoleculeInfo minfo(get_molecule_info(species));

        const auto target_space = voxel.space.lock();
        for (const auto &space : spaces_)
        {
            if (space->has_voxel(pid))
            {
                if (space != target_space)
                {
                    space->remove_voxel(pid);
                }
                return target_space->update_voxel(pid, species,
                                                  voxel.coordinate);
            }
        }

        if (!target_space->has_species(species))
        {
            target_space->make_molecular_type(species, minfo.loc);
        }
        return target_space->update_voxel(pid, species, voxel.coordinate);
    }

    bool remove_voxel(const ParticleID &pid)
    {
        for (const auto &space : spaces_)
        {
            if (space->has_voxel(pid))
                return space->remove_voxel(pid);
        }
        return false;
    }

    // Deprecated
    bool can_move(const Voxel &src, const Voxel &dst) const
    {
        if (!(src.space < dst.space) && !(dst.space < src.space))
            return src.space.lock()->can_move(src.coordinate, dst.coordinate);

        return false;
    }

    // Deprecated
    bool move(const Voxel &src, const Voxel &dst,
              const std::size_t candidate = 0)
    {
        if (!(src.space < dst.space) && !(dst.space < src.space))
            return src.space.lock()->move(src.coordinate, dst.coordinate,
                                          candidate);

        return false;
    }

    const Integer size() const { return size_; }

    const Integer3 shape() const { return get_root()->shape(); }

    /*
     * SpatiocyteWorld API
     */

    const MoleculeInfo get_molecule_info(const Species &sp) const
    {
        const auto itr = molecule_info_cache_.find(sp);
        if (itr != molecule_info_cache_.end())
        {
            return itr->second;
        }
        throw NotFound("MoleculeInfo not found");
    }

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    const MoleculeInfo get_molecule_info(const Species &sp)
    {
        boost::optional<Real> diffusion_coef;
        boost::optional<Real> radius;
        boost::optional<std::string> location;
        Shape::dimension_kind dimension =
            Shape::THREE; // The default dimension is 3

        const auto itr = molecule_info_cache_.find(sp);
        if (itr != molecule_info_cache_.end())
        {
            // return itr->second;
            // TODO: the below code is only for warning.
            //       In the future, the value should be returned immediately.
            diffusion_coef = itr->second.D;
            radius = itr->second.radius;
            location = itr->second.loc;
            dimension = itr->second.dimension;
        }
        else if (const auto model = lock_model())
        {
            const auto species_from_model(model->apply_species_attributes(sp));

            if (species_from_model.has_attribute("D"))
            {
                diffusion_coef = species_from_model.get_attribute_as<Real>("D");
            }

            if (species_from_model.has_attribute("radius"))
            {
                radius = species_from_model.get_attribute_as<Real>("radius");
            }

            if (species_from_model.has_attribute("location"))
            {
                location = species_from_model.get_attribute_as<std::string>(
                    "location");
            }

            dimension = extras::get_dimension_from_model(sp, model);
        }

        if (sp.has_attribute("D"))
        {
            const auto new_value = sp.get_attribute_as<Real>("D");
            if (diffusion_coef && *diffusion_coef != new_value)
            {
                warning("D");
                diffusion_coef = new_value;
            }
        }

        if (sp.has_attribute("radius"))
        {
            const auto new_value = sp.get_attribute_as<Real>("radius");
            if (radius && *radius != new_value)
            {
                warning("radius");
                radius = new_value;
            }
        }

        if (sp.has_attribute("location"))
        {
            const auto new_value = sp.get_attribute_as<std::string>("location");
            if (location && *location != new_value)
            {
                warning("location");
                location = new_value;
            }
        }

        if (diffusion_coef == boost::none)
        {
            diffusion_coef = 0.0;
        }

        if (radius == boost::none)
        {
            radius = voxel_radius();
        }

        if (location == boost::none)
        {
            location = "";
        }

        const MoleculeInfo info = {*radius, *diffusion_coef, *location,
                                   dimension};
        molecule_info_cache_.insert(
            molecule_info_cache_t::value_type(sp, info));
        return info;
    }

    static inline void warning(const std::string attribute)
    {
        std::cerr << "Warning: A given species has an attribute \"" << attribute
                  << "\"";
        std::cerr << ", but its value differs from that of the bound Model or "
                     "the value previously given."
                  << std::endl;
        std::cerr << "         Giving the different value from a species "
                     "attribute are deprecated."
                  << std::endl;
        std::cerr << "         An attribute of a given species will be ignored "
                     "in the future."
                  << std::endl;
    }

    // bool has_species_exact(const Species &sp) const
    // {
    //     return get_root()->has_species_exact(sp);
    // }

    void set_value(const Species &sp, const Real value);

    /**
     * create and add a new particle
     * @param p a particle
     * @return a pair of a pair of pid (a particle id) and p (a particle)
     * and bool (if it's succeeded or not)
     */
    boost::optional<ParticleID> new_particle(const Particle &p)
    {
        // ParticleID pid(sidgen_());
        // const bool is_succeeded(update_particle(pid, p));
        // return std::make_pair(get_particle(pid), is_succeeded);
        const MoleculeInfo minfo(get_molecule_info(p.species()));
        const Voxel voxel(get_voxel_nearby(p.position()));

        if (voxel.get_voxel_pool()->species().serial() != minfo.loc)
            return boost::none;

        if (boost::optional<ParticleID> pid = new_particle(p.species(), voxel))
            return *pid;

        return boost::none;
    }

    boost::optional<ParticleID> new_particle(const Species &sp,
                                             const Real3 &pos)
    {
        const MoleculeInfo info(get_molecule_info(sp));
        return new_particle(Particle(sp, pos, info.radius, info.D));
    }

    bool remove_particle(const ParticleID &pid) { return remove_voxel(pid); }

    bool update_particle(const ParticleID &pid, const Particle &p)
    {
        const MoleculeInfo minfo(get_molecule_info(p.species()));
        return update_voxel(pid, p.species(), get_voxel_nearby(p.position()));
    }

    std::vector<Species> list_species() const
    {
        std::vector<Species> list;
        for (const auto &space : spaces_)
        {
            std::vector<Species> species(space->list_species());
            list.insert(list.end(), species.begin(), species.end());
        }
        return list;
    }

    std::vector<Species> list_non_structure_species() const;
    std::vector<Species> list_structure_species() const;
    // std::vector<coordinate_type> list_coords(const Species& sp) const;

    boost::optional<ParticleID> new_particle(const Species &sp,
                                             const Voxel &voxel)
    {
        boost::shared_ptr<VoxelSpaceBase> space(voxel.space.lock());
        if (!space->has_species(sp))
        {
            const MoleculeInfo minfo(get_molecule_info(sp));
            space->make_molecular_type(sp, minfo.loc);
        }

        ParticleID pid(sidgen_());

        if (space->add_voxel(sp, pid, voxel.coordinate))
            return pid;

        return boost::none;
    }

    boost::optional<ParticleID> new_voxel_structure(const Species &sp,
                                                    const Voxel &voxel)
    {
        boost::shared_ptr<VoxelSpaceBase> space(voxel.space.lock());
        if (!space->has_species(sp))
        {
            const MoleculeInfo minfo(get_molecule_info(sp));
            space->make_structure_type(sp, minfo.loc);
        }

        ParticleID pid;

        if (space->add_voxel(sp, pid, voxel.coordinate))
            return pid;

        return boost::none;
    }

    bool add_molecules(const Species &sp, const Integer &num);
    bool add_molecules(const Species &sp, const Integer &num,
                       const boost::shared_ptr<const Shape> shape);
    Integer add_structure(const Species &sp,
                          const boost::shared_ptr<const Shape> shape);

    void remove_molecules(const Species &sp, const Integer &num);
    // void remove_molecules_exact(const Species& sp, const Integer& num);

    boost::optional<Voxel> check_neighbor(const Voxel &voxel,
                                          const std::string &loc);

    const Integer num_neighbors(const Voxel &voxel) const
    {
        return voxel.space.lock()->num_neighbors(voxel.coordinate);
    }

    const Voxel get_neighbor(const Voxel &voxel, Integer nrand) const
    {
        return Voxel(voxel.space,
                     voxel.space.lock()->get_neighbor(voxel.coordinate, nrand));
    }

    const Voxel get_neighbor_randomly(const Voxel &voxel) const;
    const Voxel get_neighbor_randomly(const Voxel &voxel,
                                      Shape::dimension_kind dimension) const;

    // bool update_molecule(coordinate_type at, Species species);

    const Species &draw_species(const Species &pttrn) const;

    boost::shared_ptr<RandomNumberGenerator> rng() const { return rng_; }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = model_.lock())
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
        const auto bound = model_.lock();
        if (!bound)
        {
            std::cerr << "Warning: Manipulating SpatiocyteWorld without "
                         "binding Model is deprecated."
                      << std::endl;
            std::cerr << "         In the future, calling this function before "
                         "binding Model will throw an exception."
                      << std::endl;
        }
        return bound;
    }

    Shape::dimension_kind get_dimension(const Species &species)
    {
        return get_molecule_info(species).dimension;
    }

    Shape::dimension_kind get_dimension(const Species &species) const
    {
        return get_molecule_info(species).dimension;
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

    static inline Integer3 calculate_shape(const Real3 &edge_lengths,
                                           const Real &voxel_radius)
    {
        return VoxelSpaceBase::calculate_shape(edge_lengths, voxel_radius,
                                               true);
    }

    static inline Real calculate_volume(const Real3 &edge_lengths,
                                        const Real &voxel_radius)
    {
        return VoxelSpaceBase::calculate_volume(edge_lengths, voxel_radius,
                                                true);
    }

protected:
    space_type get_root() const { return spaces_.at(0); }

    space_type get_space(coordinate_type &coordinate)
    {
        for (const auto &space : spaces_)
        {
            if (coordinate < space->size())
                return space;

            coordinate -= space->size();
        }
        return space_type();
    }

    Integer add_structure2(const Species &sp, const std::string &location,
                           const boost::shared_ptr<const Shape> shape);
    Integer add_structure3(const Species &sp, const std::string &location,
                           const boost::shared_ptr<const Shape> shape);
    bool is_surface_voxel(const Voxel &voxel,
                          const boost::shared_ptr<const Shape> shape) const;

    Particle gen_particle_from(const space_type &space,
                               const VoxelView &view) const
    {
        const auto species = view.species;

        const auto position = space->coordinate2position(view.voxel);
        const auto minfo_iter = molecule_info_cache_.find(species);
        if (minfo_iter != molecule_info_cache_.end())
        {
            const auto &minfo = minfo_iter->second;
            return Particle(species, position, minfo.radius, minfo.D,
                            minfo.loc);
        }
        else
        {
            return Particle(species, position, 0.0, 0.0, "");
        }
    }

private:
    template <typename T, typename ListFn, typename Fn>
    std::vector<T> map_voxels(ListFn list_f, Fn f) const
    {
        std::vector<T> list;
        for (const auto &space : spaces_)
        {
            const auto voxels(list_f(space));
            list.reserve(list.size() + voxels.size());
            for (const auto &item : voxels)
            {
                list.push_back(f(space, item));
            }
        }
        return list;
    }

    template <typename ListFn>
    std::vector<std::pair<ParticleID, Particle>>
    list_particles_private(ListFn list_fn) const
    {
        return map_voxels<std::pair<ParticleID, Particle>>(
            list_fn, [this](const space_type &space, const VoxelView &view) {
                return std::make_pair(view.pid, gen_particle_from(space, view));
            });
    }

    template <typename ListFn>
    std::vector<ParticleBase<Voxel>> list_voxels_private(ListFn list_fn) const
    {
        return map_voxels<ParticleBase<Voxel>>(
            list_fn, [this](const space_type &space, const VoxelView &view) {
                return ParticleBase<Voxel>(view.pid, view.species,
                                           Voxel(space, view.voxel));
            });
    }

protected:
    typedef utils::get_mapper_mf<Species, MoleculeInfo>::type
        molecule_info_cache_t;

    std::size_t size_;
    space_container_type spaces_;

    OneToManyMap<Voxel> interfaces_;
    OneToManyMap<Voxel> neighbors_;

    boost::shared_ptr<RandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> sidgen_;

    boost::weak_ptr<Model> model_;
    molecule_info_cache_t molecule_info_cache_;
};

inline SpatiocyteWorld *create_spatiocyte_world_cell_list_impl(
    const Real3 &edge_lengths, const Real &voxel_radius,
    const Integer3 &matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator> &rng)
{
    return new SpatiocyteWorld(
        new LatticeSpaceCellListImpl(edge_lengths, voxel_radius, matrix_sizes),
        rng);
}

inline SpatiocyteWorld *create_spatiocyte_world_vector_impl(
    const Real3 &edge_lengths, const Real &voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator> &rng)
{
    return new SpatiocyteWorld(
        new LatticeSpaceVectorImpl(edge_lengths, voxel_radius), rng);
}

inline SpatiocyteWorld *allocate_spatiocyte_world_square_offlattice_impl(
    const Real edge_length, const Species &species, const Real &voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator> &rng)
{
    OffLatticeSpace::position_container positions;
    OffLatticeSpace::coordinate_pair_list_type adjoining_pairs;

    // const std::size_t num_row(int((edge_length - (2 + sqrt(3)) *
    // voxel_radius) /
    //                          (2 * sqrt(3) * voxel_radius)) + 1);
    const std::size_t num_row(int(edge_length / (2 * sqrt(3) * voxel_radius)));
    const std::size_t num_col(int(edge_length / (2 * voxel_radius)));

    for (std::size_t row(0); row < num_row; ++row)
        for (std::size_t col(0); col < num_col; ++col)
        {
            // 2 * (row * num_col + col)
            positions.push_back(Real3(2 * col, 2 * row * sqrt(3), 0) *
                                voxel_radius);

            // 2 * (row * num_col + col + 1)
            positions.push_back(Real3(2 * col, (2 * row + 1) * sqrt(3), 0) *
                                voxel_radius);

            const int index(2 * (row * num_col + col));

            const std::size_t next_col((col + 1) % num_col);
            const std::size_t next_row((row + 1) % num_row);

            const int right(2 * (row * num_col + next_col));
            const int bottom(2 * (next_row * num_col + col));
            const int right_bottom(2 * (next_row * num_col + next_col));

            adjoining_pairs.push_back(std::make_pair(index, index + 1));
            adjoining_pairs.push_back(std::make_pair(index, right));
            adjoining_pairs.push_back(std::make_pair(index + 1, right));
            adjoining_pairs.push_back(std::make_pair(index + 1, right + 1));
            adjoining_pairs.push_back(std::make_pair(index + 1, bottom));
            adjoining_pairs.push_back(std::make_pair(index + 1, right_bottom));
        }

    OffLatticeSpace *space =
        new OffLatticeSpace(voxel_radius, species, positions, adjoining_pairs);
    space->set_lengths(Real3(2 * num_col, 2 * sqrt(3) * num_row, 2) *
                       voxel_radius);

    return new SpatiocyteWorld(space, rng);
}

} // namespace spatiocyte

} // namespace ecell4

#endif /* ECELL4_LATTICE_LATTICE_WORLD_HPP */
