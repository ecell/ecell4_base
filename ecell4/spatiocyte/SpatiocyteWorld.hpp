#ifndef ECELL4_LATTICE_LATTICE_WORLD_HPP
#define ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <sstream>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

// #include <ecell4/core/LatticeSpace.hpp>
#include <ecell4/core/LatticeSpaceCellListImpl.hpp>
#include <ecell4/core/LatticeSpaceVectorImpl.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/extras.hpp>

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
    : public Space
{
public:

    // typedef LatticeSpaceCellListImpl default_space_type;
    typedef LatticeSpaceVectorImpl default_space_type;

    typedef MoleculeInfo molecule_info_type;

    typedef LatticeSpace::coordinate_id_pair_type coordinate_id_pair_type;
    typedef LatticeSpace::coordinate_type coordinate_type;

public:

    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
        : space_(new default_space_type(edge_lengths, voxel_radius)), rng_(rng)
    {
        ; // do nothing
    }

    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius)
        : space_(new default_space_type(edge_lengths, voxel_radius))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const Real3& edge_lengths = Real3(1, 1, 1))
        : space_(new default_space_type(edge_lengths, edge_lengths[0] / 100)) //XXX: sloppy default
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const std::string filename)
        : space_(new default_space_type(Real3(1, 1, 1), 1 / 100)) //XXX: sloppy default
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }

    SpatiocyteWorld(LatticeSpace* space,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
        : space_(space), rng_(rng)
    {
        ; // do nothing
    }

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

    const Real t() const;
    void set_t(const Real& t);

    const Real3& edge_lengths() const;
    const Real volume() const;
    Integer num_species() const;
    bool has_species(const Species &sp) const;
    // bool has_species_exact(const Species &sp) const;

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;
    Integer num_voxels() const;
    Integer num_voxels(const Species& sp) const;
    Integer num_voxels_exact(const Species& sp) const;

    void set_value(const Species& sp, const Real value);

    Real get_value(const Species& sp) const
    {
        return (*space_).get_value(sp);
    }

    Real get_value_exact(const Species& sp) const
    {
        return (*space_).get_value_exact(sp);
    }

    /**
     * create and add a new particle
     * @param p a particle
     * @return a pair of a pair of pid (a particle id) and p (a particle)
     * and bool (if it's succeeded or not)
     */
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p)
    {
        // ParticleID pid(sidgen_());
        // const bool is_succeeded(update_particle(pid, p));
        // return std::make_pair(get_particle(pid), is_succeeded);
        const molecule_info_type minfo(get_molecule_info(p.species()));
        const Voxel v(
            p.species(), position2coordinate(p.position()), p.radius(), p.D(), minfo.loc);
        if ((*space_).on_structure(v))
        {
            return std::make_pair(std::make_pair(ParticleID(), p), false);
        }
        const std::pair<std::pair<ParticleID, Voxel>, bool> retval = new_voxel(v);
        return std::make_pair(std::make_pair(retval.first.first, p), retval.second);
    }

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Species& sp, const Real3& pos)
    {
        const MoleculeInfo info(get_molecule_info(sp));
        return new_particle(Particle(sp, pos, info.radius, info.D));
    }

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {
        return (*space_).get_particle(pid);
    }

    std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const
    {
        return (*space_).get_voxel(pid);
    }

    std::pair<ParticleID, Voxel> get_voxel_at(const coordinate_type& coord) const
    {
        return (*space_).get_voxel_at(coord);
    }

    bool remove_particle(const ParticleID& pid)
    {
        return (*space_).remove_particle(pid);
    }

    bool remove_voxel(const ParticleID& pid)
    {
        return (*space_).remove_voxel(pid);
    }

    bool has_voxel(const ParticleID& pid) const;
    bool has_particle(const ParticleID& pid) const;

    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles_exact(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> > list_structure_particles() const;
    std::vector<std::pair<ParticleID, Particle> > list_non_structure_particles() const;

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        const molecule_info_type minfo(get_molecule_info(p.species()));
        return update_voxel(pid, Voxel(p.species(),
            position2coordinate(p.position()), p.radius(), p.D(), minfo.loc));
    }

    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels() const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const;

    std::vector<Species> list_species() const;
    std::vector<Species> list_non_structure_species() const;
    std::vector<Species> list_structure_species() const;
    // std::vector<coordinate_type> list_coords(const Species& sp) const;

    bool has_molecule_pool(const Species& sp) const
    {
        return (*space_).has_molecule_pool(sp);
    }

    MoleculePool* find_molecule_pool(const Species& species)
    {
        return (*space_).find_molecule_pool(species);
    }

    const MoleculePool* find_molecule_pool(const Species& species) const
    {
        return (*space_).find_molecule_pool(species);
    }

    VoxelPool* find_voxel_pool(const Species& species);
    const VoxelPool* find_voxel_pool(const Species& species) const;
    VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const;

    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel(const Voxel& v);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel(const Species& sp, const coordinate_type& coord);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel_structure(const Species& sp, const coordinate_type& coord);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel_structure(const Voxel& v);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel_interface(const Species& sp, const coordinate_type& coord);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel_interface(const Voxel& v);

    bool add_molecules(const Species& sp, const Integer& num);
    bool add_molecules(const Species& sp, const Integer& num, const boost::shared_ptr<const Shape> shape);
    Integer add_structure(const Species& sp, const boost::shared_ptr<const Shape> shape);
    Integer add_interface(const Species& sp);
    Integer add_neighbors(const Species& sp, const coordinate_type center); // TODO

    void remove_molecules(const Species& sp, const Integer& num);
    // void remove_molecules_exact(const Species& sp, const Integer& num);
    bool remove_voxel(const coordinate_type coord);

    bool move(const coordinate_type& src, const coordinate_type& dest,
              const std::size_t candidate=0);
    bool can_move(const coordinate_type& src, const coordinate_type& dest) const;

    // std::pair<coordinate_type, bool> move_to_neighbor(
    //     coordinate_type coord, Integer nrand);
    // std::pair<coordinate_type, bool> move_to_neighbor(
    //     coordinate_id_pair_type& info, Integer nrand);
    // std::pair<std::pair<coordinate_id_pair_type, coordinate_type>, bool>
    //     move_to_neighbor(VoxelPool* mtype, Integer index);
    std::pair<coordinate_type, bool> move_to_neighbor(
        VoxelPool* const& from_mt, VoxelPool* const& loc,
        coordinate_id_pair_type& info, const Integer nrand);

    coordinate_type get_neighbor(coordinate_type coord, Integer nrand) const
    {
        return (*space_).get_neighbor(coord, nrand);
    }

    coordinate_type get_neighbor_boundary(
            coordinate_type coord, Integer nrand) const
    {
        return (*space_).get_neighbor_boundary(coord, nrand);
    }

    std::pair<coordinate_type, bool> check_neighbor(
            const coordinate_type coord, const std::string& loc);
    // bool update_molecule(coordinate_type at, Species species);

    const Species& draw_species(const Species& pttrn) const;

    // std::pair<std::pair<ParticleID, Voxel>, bool> place_voxel(
    //     const Species& sp, const coordinate_type& coord)
    // {
    //     const molecule_info_type& info(get_molecule_info(sp));
    //     return new_voxel(ecell4::Voxel(sp, coord, info.radius, info.D));
    // }

    // void update_voxel(const Voxel& v)
    // {
    //     (*space_).update_voxel(v);
    // }

    bool update_voxel(const ParticleID& pid, const Voxel& v)
    {
        return (*space_).update_voxel(pid, v);
    }

    Real voxel_radius() const
    {
        return (*space_).voxel_radius();
    }

    Real voxel_volume() const
    {
        return (*space_).voxel_volume();
    }

    Real unit_area() const
    {
        return (*space_).unit_area();
    }

    Real get_volume(const Species& sp) const
    {
        if (!has_species(sp) || !find_molecule_pool(sp)->is_structure())
        {
            return 0.0;
        }
        return (*space_).get_volume(sp);
    }

    Real3 actual_lengths() const
    {
        return (*space_).actual_lengths();
    }

    boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    const Integer size() const
    {
        return (*space_).size();
    }

    const Integer3 shape() const
    {
        return (*space_).shape();
    }

    const Integer inner_size() const
    {
        return (*space_).inner_size();
    }

    // TODO
    // const Integer3 inner_shape() const
    // {
    //     return (*space_).inner_shape();
    // }

    const coordinate_type inner2coordinate(const coordinate_type inner)
    {
        return (*space_).inner2coordinate(inner);
    }

    coordinate_type position2coordinate(const Real3& pos) const
    {
        return (*space_).position2coordinate(pos);
    }

    const Real3 coordinate2position(const coordinate_type& coord) const
    {
        return (*space_).coordinate2position(coord);
    }

    /**
     * temp
     */

    const molecule_info_type get_molecule_info(const VoxelPool* mt) const
    {
        const std::string loc(
            mt->location()->is_vacant() ? "" : mt->location()->species().serial());
        molecule_info_type info = {mt->radius(), mt->D(), loc};
        return info;
    }

    std::pair<ParticleID, Voxel> make_pid_voxel_pair(
        const VoxelPool* mt, const coordinate_type& coord) const
    {
        const ParticleID pid(mt->get_particle_id(coord));
        const coordinate_id_pair_type info(pid, coord);
        return make_pid_voxel_pair(mt, info);
    }

    std::pair<ParticleID, Voxel> make_pid_voxel_pair(
        const VoxelPool* mt, const coordinate_id_pair_type& info) const
    {
        const std::string loc(
            mt->location()->is_vacant() ? "" : mt->location()->species().serial());
        return std::make_pair<ParticleID, Voxel>(
            ParticleID(info.pid),
            Voxel(mt->species(), info.coordinate, mt->radius(), mt->D(), loc));
    }

    std::pair<ParticleID, Voxel> choice(const Species& sp)
    {
        const MoleculePool* mt(find_molecule_pool(sp));
        const Integer i(rng_->uniform_int(0, mt->size() - 1));
        const coordinate_id_pair_type& info(mt->at(i));
        return make_pid_voxel_pair(mt, info);
    }

    // bool on_structure(const Species& sp, const coordinate_type& coord)
    // {
    //     const molecule_info_type minfo(get_molecule_info(sp));
    //     return on_structure(
    //         Voxel(sp, coord, minfo.radius, minfo.D, minfo.loc));
    // }

    bool on_structure(const Voxel& v)
    {
        return (*space_).on_structure(v);
    }

    /*
     * HDF5 Save
     */
    void save(const std::string& filename) const
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        sidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("LatticeSpace")));
        (*space_).save_hdf5(group.get());
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
        (*space_).load_hdf5(group);
        sidgen_.load(*fin);
        rng_->load(*fin);
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
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

    /**
     * static members
     */
    static inline Real calculate_voxel_volume(const Real r)
    {
        return LatticeSpace::calculate_voxel_volume(r);
    }

    static inline Real3 calculate_hcp_lengths(const Real voxel_radius)
    {
        return LatticeSpace::calculate_hcp_lengths(voxel_radius);
    }

    static inline Integer3 calculate_shape(const Real3& edge_lengths, const Real& voxel_radius)
    {
        return LatticeSpace::calculate_shape(edge_lengths, voxel_radius, true);
    }

    static inline Real calculate_volume(const Real3& edge_lengths, const Real& voxel_radius)
    {
        return LatticeSpace::calculate_volume(edge_lengths, voxel_radius, true);
    }

protected:

    Integer add_structure2(const Species& sp, const boost::shared_ptr<const Shape> shape);
    Integer add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape);
    bool is_surface_voxel(const coordinate_type coord,
            const boost::shared_ptr<const Shape> shape) const;

protected:

    boost::scoped_ptr<LatticeSpace> space_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> sidgen_;

    boost::weak_ptr<Model> model_;
};

SpatiocyteWorld* create_spatiocyte_world_cell_list_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const Integer3& matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator>& rng);
SpatiocyteWorld* create_spatiocyte_world_vector_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng);

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
