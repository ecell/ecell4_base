#ifndef __ECELL4_LATTICE_LATTICE_WORLD_HPP
#define __ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/MolecularTypeBase.hpp>
#include <ecell4/core/LatticeSpace.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

namespace ecell4
{

namespace lattice
{

class LatticeWorld
{
public:

    typedef LatticeSpace::coordinate_type coordinate_type;
    typedef LatticeSpace::private_coordinate_type private_coordinate_type;
    typedef LatticeSpace::particle_info particle_info;

    LatticeWorld(const Position3& edge_lengths, const Real& voxel_radius,
            boost::shared_ptr<GSLRandomNumberGenerator> rng)
        : space_(edge_lengths, voxel_radius), rng_(rng)
    {
        ; // do nothing
    }

    LatticeWorld(const Position3& edge_lengths, const Real& voxel_radius)
        : space_(edge_lengths, voxel_radius)
    {
        rng_ = boost::shared_ptr<GSLRandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    LatticeWorld(const Position3& edge_lengths)
        : space_(edge_lengths, edge_lengths[0] / 100) //XXX: sloppy default
    {
        rng_ = boost::shared_ptr<GSLRandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    const Real& t() const;
    void set_t(const Real& t);

    const Position3& edge_lengths() const;
    Integer num_species() const;
    bool has_species(const Species &sp) const;

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules() const;
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;

    bool has_particle(const ParticleID& pid) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;

    bool update_particle(const ParticleID& pid, const Particle& p);

    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const;

    std::vector<Species> list_species() const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    MolecularTypeBase* get_molecular_type(const Species& species);
    MolecularTypeBase* get_molecular_type(const private_coordinate_type& coord);
    bool add_species(const Species& sp);
    std::pair<ParticleID, bool> add_molecule(const Species& sp, coordinate_type coord);
    bool add_molecules(const Species& sp, const Integer& num);
    bool remove_molecule(const coordinate_type coord);
    bool move(coordinate_type from, coordinate_type to);
    std::pair<coordinate_type, bool> move_to_neighbor(coordinate_type coord, Integer nrand);
    std::pair<coordinate_type, bool> move_to_neighbor(particle_info& info, Integer nrand);
    bool update_molecule(coordinate_type at, Species species);

    Real voxel_radius() const
    {
        return space_.voxel_radius();
    }

    boost::shared_ptr<GSLRandomNumberGenerator> rng()
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

    coordinate_type global2coord(const Global& global) const;
    const Global coord2global(coordinate_type coord) const;

    coordinate_type private2coord(const private_coordinate_type&
            private_coord) const;
    private_coordinate_type coord2private(const coordinate_type&
            coord) const;

    /*
     * HDF5 Save
     */
    void save(const std::string& filename) const
    {
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename, H5F_ACC_TRUNC));
        rng_->save(fout.get());
        sidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("LatticeSpace")));
        space_.save(group.get());
    }

    void load(const std::string& filename)
    {
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename, H5F_ACC_RDONLY));
        const H5::Group group(fin->openGroup("LatticeSpace"));
        space_.load(group);
        sidgen_.load(*fin);
        rng_->load(*fin);
    }

protected:

    LatticeSpace space_;
    boost::shared_ptr<GSLRandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> sidgen_;

};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_WORLD_HPP */
