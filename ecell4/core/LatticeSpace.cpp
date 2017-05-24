#include "LatticeSpace.hpp"
// #include <cmath>
// #include <sstream>
// #include <algorithm>

namespace ecell4
{

bool LatticeSpace::can_move(const coordinate_type& src,
        const coordinate_type& dest) const
{
    throw NotImplemented("can_move is not implemented.");
    return false;
}

bool LatticeSpace::make_structure_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    throw NotImplemented("make_structure_type is not implemented.");
    return false;
}

bool LatticeSpace::make_interface_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    throw NotImplemented("make_interface_type is not implemented.");
    return false;
}

std::vector<std::pair<ParticleID, Particle> >
LatticeSpace::list_particles() const
{
    const std::vector<std::pair<ParticleID, Voxel> > voxels(list_voxels());

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(voxels.size());
    for (std::vector<std::pair<ParticleID, Voxel> >::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate()));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
LatticeSpace::list_particles(const Species& sp) const
{
    const std::vector<std::pair<ParticleID, Voxel> > voxels(list_voxels(sp));

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(voxels.size());
    for (std::vector<std::pair<ParticleID, Voxel> >::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate()));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
LatticeSpace::list_particles_exact(const Species& sp) const
{
    const std::vector<std::pair<ParticleID, Voxel> >
        voxels(list_voxels_exact(sp));

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(voxels.size());
    for (std::vector<std::pair<ParticleID, Voxel> >::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate()));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

std::pair<ParticleID, Particle> LatticeSpace::get_particle(const ParticleID& pid) const
{
    const Voxel v(get_voxel(pid).second);
    return std::make_pair(pid, Particle(
        v.species(), coordinate2position(v.coordinate()), v.radius(), v.D()));
}

} // ecell4
