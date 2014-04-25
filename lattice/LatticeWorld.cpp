#include "LatticeWorld.hpp"

namespace ecell4
{

namespace lattice
{

const Real& LatticeWorld::t() const
{
    return space_.t();
}

void LatticeWorld::set_t(const Real& t)
{
    space_.set_t(t);
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

std::vector<std::pair<ParticleID, Voxel> >
    LatticeWorld::list_voxels(const Species& sp) const
{
    return space_.list_voxels(sp);
}

std::vector<LatticeWorld::coordinate_type> LatticeWorld::list_coords(const Species& sp) const
{
    return space_.list_coords(sp);
}

MolecularTypeBase* LatticeWorld::find_molecular_type(const Species& species)
{
    return space_.find_molecular_type(species);
}

MolecularTypeBase* LatticeWorld::get_molecular_type(const private_coordinate_type& coord)
{
    return space_.get_molecular_type(coord);
}

// bool LatticeWorld::register_species(const Species& sp)
// {
//     return space_.register_species(sp);
// }

// std::pair<ParticleID, bool> LatticeWorld::add_molecule(const Species& sp, const private_coordinate_type& coord)
// {
//     ParticleID pid(sidgen_());
//     return std::pair<ParticleID, bool>(pid, space_.add_molecule(sp, coord, pid));
// }

std::pair<ParticleID, bool> LatticeWorld::new_voxel_private(const Voxel& v)
{
    ParticleID pid(sidgen_());
    return std::make_pair(pid, space_.update_voxel_private(pid, v));
}

bool LatticeWorld::add_molecules(const Species& sp, const Integer& num)
{
    const LatticeWorld::molecule_info_type info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const coordinate_type coord(rng()->uniform_int(0, space_.size() - 1));
        if (new_voxel_private(
            Voxel(sp, coord2private(coord), info.radius, info.D)).second)
        {
            ++count;
        }
    }
    return true;
}

bool LatticeWorld::remove_molecule(const coordinate_type coord)
{
    return space_.remove_molecule(coord);
}

bool LatticeWorld::move(coordinate_type from, coordinate_type to)
{
    return space_.move(from, to);
}

std::pair<LatticeWorld::coordinate_type, bool> LatticeWorld::move_to_neighbor(
        coordinate_type coord, Integer nrand)
{
    return space_.move_to_neighbor(coord, nrand);
}

std::pair<LatticeWorld::coordinate_type, bool> LatticeWorld::move_to_neighbor(
        particle_info& info, Integer nrand)
{
    return space_.move_to_neighbor(info, nrand);
}

std::pair<std::pair<LatticeWorld::particle_info,
    LatticeWorld::private_coordinate_type>, bool>
LatticeWorld::move_to_neighbor(MolecularTypeBase* mtype, Integer index)
{
    const Integer rnd(rng()->uniform_int(0,11));
    particle_info& info(mtype->at(index));
    std::pair<private_coordinate_type, bool> neighbor(
            space_.move_to_neighbor(info, rnd));
    return std::make_pair(std::make_pair(info, neighbor.first), neighbor.second);
}

std::pair<LatticeWorld::private_coordinate_type, bool> LatticeWorld::check_neighbor(
        const private_coordinate_type coord)
{
    const Integer rnd(rng()->uniform_int(0,11));
    const private_coordinate_type neighbor(space_.get_neighbor(coord, rnd));
    bool flg = get_molecular_type(neighbor)->is_vacant();
    return std::make_pair(neighbor, flg);
}

// bool LatticeWorld::update_molecule(coordinate_type at, Species species)
// {
//     return space_.update_molecule(at, species);
// }

LatticeWorld::coordinate_type LatticeWorld::global2coord(const Global& global) const
{
    return space_.global2coord(global);
}

const Global LatticeWorld::coord2global(coordinate_type coord) const
{
    return space_.coord2global(coord);
}

LatticeWorld::coordinate_type LatticeWorld::private2coord(
        const private_coordinate_type& private_coord) const
{
    return space_.private2coord(private_coord);
}

LatticeWorld::private_coordinate_type LatticeWorld::coord2private(
        const coordinate_type& coord) const
{
    return space_.coord2private(coord);
}

} // lattice

} // ecell4
