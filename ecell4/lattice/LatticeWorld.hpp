#ifndef __ECELL4_LATTICE_LATTICE_WORLD_HPP
#define __ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/LatticeSpace.hpp>
#include <ecell4/core/MolecularTypeBase.hpp>
#include <ecell4/core/MolecularType.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Shape.hpp>

namespace ecell4
{

namespace lattice
{

struct MoleculeInfo
{
    const Real radius;
    const Real D;
    const std::string loc;
};

class LatticeWorld
    : public Space
{
public:

    typedef LatticeSpace::coordinate_type coordinate_type;
    typedef LatticeSpace::private_coordinate_type private_coordinate_type;
    typedef LatticeSpace::particle_info particle_info;
    typedef LatticeSpace::spmap spmap;
    typedef MoleculeInfo molecule_info_type;

    typedef std::map<std::string, Shape::dimension_kind> dimension_map_type;

public:

    LatticeWorld(const Real3& edge_lengths, const Real& voxel_radius,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
        : space_(edge_lengths, voxel_radius), rng_(rng)
    {
        ; // do nothing
    }

    LatticeWorld(const Real3& edge_lengths, const Real& voxel_radius)
        : space_(edge_lengths, voxel_radius)
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    LatticeWorld(const Real3& edge_lengths = Real3(1, 1, 1))
        : space_(edge_lengths, edge_lengths[0] / 100) //XXX: sloppy default
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    LatticeWorld(const std::string filename)
        : space_(Real3(1, 1, 1), 1 / 100) //XXX: sloppy default
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
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
    Integer num_molecules() const;
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;
    Integer num_voxels() const;
    Integer num_voxels(const Species& sp) const;
    Integer num_voxels_exact(const Species& sp) const;

    Real get_value(const Species& sp) const
    {
        return space_.get_value(sp);
    }

    Real get_value_exact(const Species& sp) const
    {
        return space_.get_value_exact(sp);
    }

    const spmap& molecular_types() const
    {
        return space_.molecular_types();
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
        ParticleID pid(sidgen_());
        // if (has_particle(pid))
        // {
        //     throw AlreadyExists("particle already exists");
        // }
        const bool is_succeeded(space_.update_particle(pid, p));
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
        return space_.get_particle(pid);
    }

    std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const
    {
        return space_.get_voxel(pid);
    }

    bool remove_particle(const ParticleID& pid)
    {
        return space_.remove_particle(pid);
    }

    bool remove_voxel(const ParticleID& pid)
    {
        return space_.remove_voxel(pid);
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
        list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const;

    std::vector<Species> list_species() const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
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
    // std::pair<coordinate_type, bool> move_to_neighbor(coordinate_type coord, Integer nrand);
    // std::pair<coordinate_type, bool> move_to_neighbor(particle_info& info, Integer nrand);
    std::pair<std::pair<particle_info, private_coordinate_type>, bool>
        move_to_neighbor(MolecularTypeBase* mtype, Integer index);

    std::pair<private_coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info& info, const Integer nrand);

    private_coordinate_type get_neighbor(
            private_coordinate_type private_coord, Integer nrand) const
    {
        return space_.get_neighbor(private_coord, nrand);
    }

    std::pair<private_coordinate_type, bool> check_neighbor_private(
            const private_coordinate_type coord);
    // bool update_molecule(coordinate_type at, Species species);

    const Species& draw_species(const Species& pttrn) const;

    std::pair<std::pair<ParticleID, Voxel>, bool> place_voxel_private(const Species& sp, const private_coordinate_type& coord)
    {
        const molecule_info_type& info(get_molecule_info(sp));
        return new_voxel_private(ecell4::Voxel(sp, coord, info.radius, info.D));
    }

    bool update_voxel(const ParticleID& pid, const Voxel& v)
    {
        return space_.update_voxel(pid, v);
    }

    bool update_voxel_private(const Voxel& v)
    {
        return space_.update_voxel_private(v);
    }

    bool update_voxel_private(const ParticleID& pid, const Voxel& v)
    {
        return space_.update_voxel_private(pid, v);
    }

    Real voxel_radius() const
    {
        return space_.voxel_radius();
    }

    Real voxel_volume() const
    {
        const Real r(voxel_radius());
        return 4.0 * sqrt(2.0) * r * r * r;
    }

    boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    const Integer col_size() const
    {
        return space_.col_size();
    }

    const Integer row_size() const
    {
        return space_.row_size();
    }

    const Integer layer_size() const
    {
        return space_.layer_size();
    }

    const Integer size() const
    {
        return space_.size();
    }

    coordinate_type position2coordinate(const Real3& pos) const
    {
        return space_.position2coordinate(pos);
    }

    const Real3 coordinate2position(const coordinate_type& coord) const
    {
        return space_.coordinate2position(coord);
    }

    const Real3 private2position(const private_coordinate_type& coord) const
    {
        return space_.coordinate2position(private2coord(coord));
    }

    const Real3 global2position(const Integer3& global) const
    {
        return space_.global2position(global);
    }

    const Integer3 position2global(const Real3& pos) const
    {
        return space_.position2global(pos);
    }

    coordinate_type global2coord(const Integer3& global) const;
    const Integer3 coord2global(coordinate_type coord) const;
    coordinate_type private2coord(const private_coordinate_type&
            private_coord) const;
    private_coordinate_type coord2private(const coordinate_type&
            coord) const;

    Shape::dimension_kind get_dimension_kind(const std::string& name) const;

    /*
     * HDF5 Save
     */
    void save(const std::string& filename) const
    {
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        sidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("LatticeSpace")));
        space_.save(group.get());
    }

    void load(const std::string& filename)
    {
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
        const H5::Group group(fin->openGroup("LatticeSpace"));
        space_.load(group);
        sidgen_.load(*fin);
        rng_->load(*fin);
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to LatticeWorld"
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

    LatticeSpace space_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> sidgen_;
    dimension_map_type dimension_map_;

    boost::weak_ptr<Model> model_;
};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_WORLD_HPP */
