#ifndef __ECELL4_LATTICE_LATTICE_WORLD_HPP
#define __ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

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

    LatticeWorld(const Position3& edge_lengths, const Real& voxel_radius,
            boost::shared_ptr<GSLRandomNumberGenerator> rng)
        : space_(edge_lengths, voxel_radius, true), t_(0), rng_(rng)
    {
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

    std::vector<Species> list_species() const;
    std::vector<Coord> list_coords(const Species& sp) const;
    MolecularTypeBase* get_molecular_type(const Species& species);
    MolecularTypeBase* get_molecular_type(Coord coord);
    bool add_species(const Species& sp);
    bool add_molecule(const Species& sp, Coord coord);
    bool add_molecules(const Species& sp, const Integer& num);
    bool remove_molecule(const Coord coord);
    std::pair<Coord, bool> move(Coord from, Coord to);
    std::pair<Coord, bool> move_to_neighbor(Coord coord, Integer nrand);
    bool update_molecule(Coord at, Species species);

    Real normalized_voxel_radius() const
    {
        return space_.normalized_voxel_radius();
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

    Coord global2coord(const Global& global) const;
    const Global coord2global(Coord coord) const;

    /*
     * HDF5 Save
     */
    void save(H5::H5File* fout, const std::string& hdf5path) const
    {
        space_.save(fout, hdf5path);
    }

protected:

    LatticeSpace space_;
    Real t_;
    boost::shared_ptr<GSLRandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> sidgen_;

};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_WORLD_HPP */
