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
    const Integer NUM_COLROW(space_.num_colrow());
    const Integer NUM_ROW(space_.num_row());
    const bool odd_col((coord % NUM_COLROW / NUM_ROW) & 1);
    const bool odd_lay((coord / NUM_COLROW) & 1);
    switch(nrand)
    {
    case 1:
        return coord+1;
    case 2:
        return coord+(odd_col^odd_lay)-NUM_ROW-1;
    case 3:
        return coord+(odd_col^odd_lay)-NUM_ROW;
    case 4:
        return coord+(odd_col^odd_lay)+NUM_ROW-1;
    case 5:
        return coord+(odd_col^odd_lay)+NUM_ROW;
    case 6:
        return coord-NUM_COLROW-(odd_col&odd_lay)-(NUM_ROW&(-!odd_lay));
    case 7:
        return coord-NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 8:
        return coord-NUM_COLROW+(odd_col&!odd_lay)+(NUM_ROW&(-odd_lay));
    case 9:
        return coord+NUM_COLROW-(odd_col&odd_lay)-(NUM_ROW&(-!odd_lay));
    case 10:
        return coord+NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 11:
        return coord+NUM_COLROW+(odd_col&!odd_lay)+(NUM_ROW&(-odd_lay));
    }
    return coord-1;
}

bool LatticeWorld::add_species(const Species& sp)
{
    return space_.add(sp);
}

bool LatticeWorld::add_molecule(const Species& sp, Coord coord)
{
    ParticleID pid(sidgen_());
    return space_.add(sp, coord, pid);
}

bool LatticeWorld::add_molecules(const Species& sp, const Integer& num)
{
    // TODO
    if (has_species(sp))
    {
        add_species(sp);
    }
    for (Integer i(0); i < num; ++i)
    {
        Coord coord(0);
        add_molecule(sp, coord);
    }
    return true;
}

bool LatticeWorld::move(Coord from, Coord to)
{
    space_.move(from, to);
}

bool LatticeWorld::react(Coord at, Species species)
{
    space_.react(at, species);
}

} // lattice

} // ecell4
