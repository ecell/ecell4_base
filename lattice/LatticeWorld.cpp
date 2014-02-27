#include "LatticeWorld.hpp"

namespace ecell4
{

namespace lattice
{

const Real& LatticeWorld::t() const
{
    return t_;
}

void LatticeWorld::set_t(const Real& t)
{
    if (t < 0.0)
    {
        throw std::invalid_argument("the time must be positive.");
    }
    t_ = t;
}

const Position3& LatticeWorld::edge_lengths() const
{
    return space_.edge_lengths();
}

Integer LatticeWorld::num_species() const
{
    return space_.num_species();
}

bool LatticeWorld::has_species(const Species &sp) const
{
    return space_.has_species(sp);
}

Integer LatticeWorld::num_molecules(const Species& sp) const
{
    return space_.num_molecules(sp);
}

Integer LatticeWorld::num_molecules() const
{
    return space_.num_particles();
}

Integer LatticeWorld::num_particles(const Species& sp) const
{
    return space_.num_particles(sp);
}

bool LatticeWorld::has_particle(const ParticleID& pid) const
{
    return space_.has_particle(pid);
}

std::vector<std::pair<ParticleID, Particle> >
LatticeWorld::list_particles() const
{
    return space_.list_particles();
}

std::vector<std::pair<ParticleID, Particle> >
LatticeWorld::list_particles(const Species& sp) const
{
    return space_.list_particles(sp);
}


bool LatticeWorld::update_particle(const ParticleID& pid, const Particle& p)
{
    return space_.update_particle(pid, p);
}


std::vector<Species> LatticeWorld::list_species() const
{
    return space_.list_species();
}

std::vector<Coord> LatticeWorld::list_coords(const Species& sp) const
{
    return space_.list_coords(sp);
}

MolecularTypeBase* LatticeWorld::get_molecular_type(const Species& species)
{
    return space_.get_molecular_type(species);
}

MolecularTypeBase* LatticeWorld::get_molecular_type(Coord coord)
{
    return space_.get_molecular_type(coord);
}

Coord LatticeWorld::get_neighbor(Coord coord, Integer nrand) const
{
	return space_.get_neighbor(coord, nrand);
}

bool LatticeWorld::add_species(const Species& sp)
{
    return space_.add_species(sp);
}

bool LatticeWorld::add_molecule(const Species& sp, Coord coord) throw(std::out_of_range)
{
    ParticleID pid(sidgen_());
    return space_.add_molecule(sp, coord, pid);
}

bool LatticeWorld::add_molecules(const Species& sp, const Integer& num)
{
    // TODO
    if (has_species(sp))
    {
        add_species(sp);
    }
    Integer count(0);
    while(count < num) {
        Integer coord(rng()->uniform_int(0,space_.size()-1));;
        if (add_molecule(sp, coord))
            ++count;
    }
    return true;
}

bool LatticeWorld::move(Coord from, Coord to) throw(std::out_of_range)
{
    space_.move(from, to);
}

bool LatticeWorld::move_to_neighbor(Coord coord, Integer nrand)
{
	space_.move_to_neighbor(coord, nrand);
}

bool LatticeWorld::react(Coord at, Species species) throw(std::out_of_range)
{
    space_.react(at, species);
}

Coord LatticeWorld::global2coord(const Global& global) const
{
    return space_.global2coord(global);
}

const Global LatticeWorld::coord2global(Coord coord) const
{
    return space_.coord2global(coord);
}

} // lattice

} // ecell4
