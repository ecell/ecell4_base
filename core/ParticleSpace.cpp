#include <cmath>
#include <stdexcept>

#include "exceptions.hpp"
#include "ParticleSpace.hpp"


namespace ecell4
{

Integer ParticleSpaceVectorImpl::num_particles() const
{
    return static_cast<Integer>(particles_.size());
}

Integer ParticleSpaceVectorImpl::num_particles(Species const& sp) const
{
    return static_cast<Integer>(list_particles(sp).size());
}

bool ParticleSpaceVectorImpl::has_particle(ParticleID const& pid) const
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    return (i != index_map_.end());
}

bool ParticleSpaceVectorImpl::update_particle(
    ParticleID const& pid, Particle const& p)
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        particle_container_type::size_type idx(particles_.size());
        index_map_[pid] = idx;
        particles_.push_back(std::make_pair(pid, p));
        return false;
    }
    else
    {
        particles_[(*i).second] = std::make_pair(pid, p);
        return true;
    }
}

void ParticleSpaceVectorImpl::remove_particle(ParticleID const& pid)
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        throw NotFound("particle not found");
    }

    particle_map_type::mapped_type
        idx((*i).second),last_idx(particles_.size() - 1);
    if (idx != last_idx)
    {
        std::pair<ParticleID, Particle> const& last(particles_[last_idx]);
        particles_[idx] = last;
        index_map_[last.first] = idx;
    }

    particles_.pop_back();
    index_map_.erase((*i).first);
}

std::pair<ParticleID, Particle> ParticleSpaceVectorImpl::get_particle(
    ParticleID const& pid) const
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        throw NotFound("particle not found");
    }

    return particles_[(*i).second];
}

std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl::list_particles() const
{
    return particles_;
}

std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl::list_particles(Species const& species) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if ((*i).second.species() == species)
        {
            retval.push_back(*i);
        }
    }

    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::list_particles_within_radius(
    Position3 const& pos, Real const& radius) const
{
    const Real rsq(gsl_pow_2(radius));
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dsq(distance_sq((*i).second.position(), pos));
        if (dsq <= rsq)
        {
            retval.push_back(std::make_pair(*i, std::sqrt(dsq)));
        }
    }

    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::list_particles_within_radius(
    Position3 const& pos, Real const& radius, ParticleID const& ignore) const
{
    const Real rsq(gsl_pow_2(radius));
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dsq(distance_sq((*i).second.position(), pos));
        if (dsq <= rsq)
        {
            if ((*i).first != ignore)
            {
                retval.push_back(std::make_pair(*i, std::sqrt(dsq)));
            }
        }
    }

    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::list_particles_within_radius(
    Position3 const& pos, Real const& radius,
    ParticleID const& ignore1, ParticleID const& ignore2) const
{
    const Real rsq(gsl_pow_2(radius));
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dsq(distance_sq((*i).second.position(), pos));
        if (dsq <= rsq)
        {
            if ((*i).first != ignore1 && (*i).first != ignore2)
            {
                retval.push_back(std::make_pair(*i, std::sqrt(dsq)));
            }
        }
    }

    return retval;
}

void ParticleSpaceVectorImpl::set_edge_lengths(Position3 const& edge_lengths)
{
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        if (edge_lengths[dim] <= 0)
        {
            throw std::invalid_argument("the edge length must be positive.");
        }
    }

    edge_lengths_ = edge_lengths;
}

} // ecell4
