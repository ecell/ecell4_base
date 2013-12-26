#ifndef __ECELL4_LATTICE_LATTICE_WORLD_HPP
#define __ECELL4_LATTICE_LATTICE_WORLD_HPP

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

    LatticeWorld(const Position3& edge_lengths,
            boost::shared_ptr<GSLRandomNumberGenerator> rng)
        : t_(0), rng_(rng), space_(edge_lengths)
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
    std::vector<std::pair<ParticleID, Particle> > list_particles(const Species& sp) const;

    bool update_particle(const ParticleID& pid, const Particle& p);

    std::vector<Species> list_species() const;
    MolecularTypeBase* get_molecular_type(const Species& species);
    MolecularTypeBase* get_molecular_type(Integer coord);
    Integer get_neighbor(Integer coord, Integer nrand) const;
    bool add_species(const Species& sp);
    bool add_molecule(const Species& sp, Coord coord);
    bool add_molecules(const Species& sp, const Integer& num);
    bool move(Coord from, Coord to);
    bool react(Coord at, Species species);

    inline boost::shared_ptr<GSLRandomNumberGenerator> rng()
    {
        return rng_;
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
