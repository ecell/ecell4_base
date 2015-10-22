#ifndef __ECELL4_LATTICE_LATTICE_WORLD_HPP
#define __ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/LatticeSpace.hpp>
#include <ecell4/core/LatticeSpaceCellListImpl.hpp>
#include <ecell4/core/MolecularTypeBase.hpp>
#include <ecell4/core/MolecularType.hpp>
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

    typedef LatticeSpace::particle_info_type particle_info_type;
    typedef LatticeSpace::coordinate_type coordinate_type;
    typedef LatticeSpace::private_coordinate_type private_coordinate_type;

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

        if (with_D && with_radius)
        {
            radius = std::atof(sp.get_attribute("radius").c_str());
            D = std::atof(sp.get_attribute("D").c_str());

            if (with_loc)
            {
                loc = sp.get_attribute("location");
            }
        }
        else
        {
            if (with_D)
            {
                D = std::atof(sp.get_attribute("D").c_str());
            }

            if (with_radius)
            {
                radius = std::atof(sp.get_attribute("radius").c_str());
            }

            if (with_loc)
            {
                loc = sp.get_attribute("location");
            }

            if (boost::shared_ptr<Model> bound_model = lock_model())
            {
                Species attributed(bound_model->apply_species_attributes(sp));
                if (!with_D && attributed.has_attribute("D"))
                {
                    D = std::atof(attributed.get_attribute("D").c_str());
                }
                if (!with_radius && attributed.has_attribute("radius"))
                {
                    radius = std::atof(
                        attributed.get_attribute("radius").c_str());
                }
                if (!with_loc && attributed.has_attribute("location"))
                {
                    loc = attributed.get_attribute("location");
                }
            }
        }

        MoleculeInfo info = {radius, D, loc};
        return info;
    }

    const Real& t() const;
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

    Real get_value(const Species& sp) const
    {
        return (*space_).get_value(sp);
    }

    Real get_value_exact(const Species& sp) const
    {
        return (*space_).get_value_exact(sp);
    }

    // const spmap& molecular_types() const
    // {
    //     return (*space_).molecular_types();
    // }

    /**
     * create and add a new particle
     * @param p a particle
     * @return a pair of a pair of pid (a particle id) and p (a particle)
     * and bool (if it's succeeded or not)
     */
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p)
    {
        ParticleID pid(sidgen_());
        // if (has_particle(pid))
        // {
        //     throw AlreadyExists("particle already exists");
        // }
        const bool is_succeeded((*space_).update_particle(pid, p));
        return std::make_pair(get_particle(pid), is_succeeded);
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

    std::pair<ParticleID, Voxel> get_voxel(const coordinate_type& coord) const
    {
        return (*space_).get_voxel(coord);
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

    bool update_particle(const ParticleID& pid, const Particle& p);

    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels() const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const;

    std::vector<Species> list_species() const;
    // std::vector<coordinate_type> list_coords(const Species& sp) const;
    MolecularTypeBase* find_molecular_type(const Species& species);
    MolecularTypeBase* get_molecular_type_private(const private_coordinate_type& coord);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel(const Voxel& v);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel(const Species& sp, const coordinate_type& coord);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel_private(const Voxel& v);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel_private(const Species& sp, const private_coordinate_type& coord);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_voxel_structure(const Voxel& v);
    bool add_molecules(const Species& sp, const Integer& num);
    bool add_molecules(const Species& sp, const Integer& num, const boost::shared_ptr<const Shape> shape);
    Integer add_structure(const Species& sp, const boost::shared_ptr<const Shape> shape);
    Integer add_neighbors(const Species& sp,
            const private_coordinate_type center); // TODO
    void remove_molecules(const Species& sp, const Integer& num);
    // void remove_molecules_exact(const Species& sp, const Integer& num);
    bool remove_voxel_private(const private_coordinate_type coord);
    bool move(coordinate_type from, coordinate_type to);
    bool move_private(const private_coordinate_type& src, const private_coordinate_type& dest);
    bool can_move(const private_coordinate_type& src, const private_coordinate_type& dest) const;
    // std::pair<coordinate_type, bool> move_to_neighbor(coordinate_type coord, Integer nrand);
    // std::pair<coordinate_type, bool> move_to_neighbor(particle_info_type& info, Integer nrand);
    // std::pair<std::pair<particle_info_type, private_coordinate_type>, bool>
    //     move_to_neighbor(MolecularTypeBase* mtype, Integer index);

    std::pair<private_coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand);

    coordinate_type get_neighbor(coordinate_type coord, Integer nrand) const
    {
        return (*space_).get_neighbor(coord, nrand);
    }

    private_coordinate_type get_neighbor_private(
            private_coordinate_type private_coord, Integer nrand) const
    {
        return (*space_).get_neighbor_private(private_coord, nrand);
    }

    std::pair<private_coordinate_type, bool> check_neighbor_private(
            const private_coordinate_type coord, const std::string& loc);
    // bool update_molecule(coordinate_type at, Species species);

    const Species& draw_species(const Species& pttrn) const;

    // std::pair<std::pair<ParticleID, Voxel>, bool> place_voxel_private(const Species& sp, const private_coordinate_type& coord)
    // {
    //     const molecule_info_type& info(get_molecule_info(sp));
    //     return new_voxel_private(ecell4::Voxel(sp, coord, info.radius, info.D));
    // }

    bool update_voxel(const ParticleID& pid, const Voxel& v)
    {
        return (*space_).update_voxel(pid, v);
    }

    void update_voxel_private(const Voxel& v)
    {
        (*space_).update_voxel_private(v);
    }

    bool update_voxel_private(const ParticleID& pid, const Voxel& v)
    {
        return (*space_).update_voxel_private(pid, v);
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

    Real get_volume() const
    {
        return (*space_).get_volume();
    }

    Real3 actual_lengths() const
    {
        return (*space_).actual_lengths();
    }

    boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    const Integer col_size() const
    {
        return (*space_).col_size();
    }

    const Integer row_size() const
    {
        return (*space_).row_size();
    }

    const Integer layer_size() const
    {
        return (*space_).layer_size();
    }

    const Integer size() const
    {
        return (*space_).size();
    }

    const Integer3 shape() const
    {
        return (*space_).shape();
    }

    coordinate_type position2coordinate(const Real3& pos) const
    {
        return (*space_).position2coordinate(pos);
    }

    const Real3 coordinate2position(const coordinate_type& coord) const
    {
        return (*space_).coordinate2position(coord);
    }

    private_coordinate_type position2private(const Real3& pos) const
    {
        return (*space_).position2private(pos);
    }

    const Real3 private2position(const private_coordinate_type& coord) const
    {
        return (*space_).coordinate2position(private2coord(coord));
    }

    coordinate_type private2coord(
        const private_coordinate_type& private_coord) const
    {
        return (*space_).private2coord(private_coord);
    }

    private_coordinate_type coord2private(
        const coordinate_type& coord) const
    {
        return (*space_).coord2private(coord);
    }

    const Real3 global2position(const Integer3& global) const
    {
        return (*space_).global2position(global);
    }

    const Integer3 position2global(const Real3& pos) const
    {
        return (*space_).position2global(pos);
    }

    coordinate_type global2coord(const Integer3& global) const
    {
        return (*space_).global2coord(global);
    }

    const Integer3 coord2global(coordinate_type coord) const
    {
        return (*space_).coord2global(coord);
    }

    private_coordinate_type global2private(const Integer3& global) const
    {
        return (*space_).global2private_coord(global);
    }

    const Integer3 private2global(private_coordinate_type coord) const
    {
        return (*space_).private_coord2global(coord);
    }

    Shape::dimension_kind get_dimension_kind(const std::string& name) const;

    /**
     * temp
     */

    const molecule_info_type get_molecule_info(const MolecularTypeBase* mt) const
    {
        const std::string loc(
            mt->location()->is_vacant() ? "" : mt->location()->species().serial());
        molecule_info_type info = {mt->radius(), mt->D(), loc};
        return info;
    }

    std::pair<ParticleID, Voxel> make_pid_voxel_pair(
        const MolecularTypeBase* mt, const private_coordinate_type& private_coord) const
    {
        const ParticleID pid(
            mt->with_voxels()
                ? mt->find_particle_id(private_coord)
                : ParticleID());
        const particle_info_type info(
            std::make_pair(private_coord, pid));
        return make_pid_voxel_pair(mt, info);
    }

    std::pair<ParticleID, Voxel> make_pid_voxel_pair(
        const MolecularTypeBase* mt, const particle_info_type& info) const
    {
        const std::string loc(
            mt->location()->is_vacant() ? "" : mt->location()->species().serial());
        const coordinate_type coord(private2coord(info.first));

        return std::make_pair<ParticleID, Voxel>(
            ParticleID(info.second()), Voxel(mt->species(), coord, mt->radius(), mt->D(), loc));
    }

    std::pair<ParticleID, Voxel> choice(const Species& sp)
    {
        MolecularTypeBase* mt(find_molecular_type(sp));
        if (!mt->with_voxels())
        {
            throw NotSupported(
                "choice for a Species with no voxel is not supporeted.");
        }
        const Integer i(rng_->uniform_int(0, mt->size() - 1));
        const particle_info_type& info(mt->at(i));
        return make_pid_voxel_pair(mt, info);
    }

    bool on_structure(const Species& sp, const coordinate_type& coord)
    {
        const private_coordinate_type private_coord(coord2private(coord));
        const molecule_info_type minfo(get_molecule_info(sp));
        return on_structure(
            Voxel(sp, private_coord, minfo.radius, minfo.D, minfo.loc));
    }

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
        (*space_).save(group.get());
        extras::save_version_information(fout.get(), "ecell4-lattice-0.0-1");
#else
        throw NotSupported("not supported yet.");
#endif
    }

    void load(const std::string& filename)
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
        const H5::Group group(fin->openGroup("LatticeSpace"));
        (*space_).load(group);
        sidgen_.load(*fin);
        rng_->load(*fin);
#else
        throw NotSupported("not supported yet.");
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

protected:

    Integer add_structure2(const Species& sp, const boost::shared_ptr<const Shape> shape);
    Integer add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape);
    bool is_surface_voxel(const Integer3& g, const boost::shared_ptr<const Shape> shape) const;

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

#endif /* __ECELL4_LATTICE_LATTICE_WORLD_HPP */
