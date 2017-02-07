#include <cmath>
#include <stdexcept>

#include "exceptions.hpp"
#include "Context.hpp"
#include "comparators.hpp"
#include "ParticleSpace.hpp"


namespace ecell4
{

Integer ParticleSpaceVectorImpl::num_particles() const
{
    return static_cast<Integer>(particles_.size());
}

Integer ParticleSpaceVectorImpl::num_particles(const Species& sp) const
{
    return static_cast<Integer>(list_particles(sp).size());
}

Integer ParticleSpaceVectorImpl::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        retval += sexp.count((*i).second.species());
    }
    return retval;
}

Integer ParticleSpaceVectorImpl::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

Integer ParticleSpaceVectorImpl::num_particles_exact(const Species& sp) const
{
    return static_cast<Integer>(list_particles_exact(sp).size());
}

bool ParticleSpaceVectorImpl::has_particle(const ParticleID& pid) const
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    return (i != index_map_.end());
}

/**
 * update or add a particle.
 * @return true if adding a new particle
 */
bool ParticleSpaceVectorImpl::update_particle(
    const ParticleID& pid, const Particle& p)
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        particle_container_type::size_type idx(particles_.size());
        index_map_[pid] = idx;
        particles_.push_back(std::make_pair(pid, p));
        return true;
    }
    else
    {
        particles_[(*i).second] = std::make_pair(pid, p);
        return false;
    }
}

void ParticleSpaceVectorImpl::remove_particle(const ParticleID& pid)
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
        const std::pair<ParticleID, Particle>& last(particles_[last_idx]);
        particles_[idx] = last;
        index_map_[last.first] = idx;
    }

    particles_.pop_back();
    index_map_.erase((*i).first);
}

std::pair<ParticleID, Particle> ParticleSpaceVectorImpl::get_particle(
    const ParticleID& pid) const
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
ParticleSpaceVectorImpl::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    SpeciesExpressionMatcher sexp(sp);

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if (sexp.match((*i).second.species()))
        {
            retval.push_back(*i);
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl::list_particles_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if ((*i).second.species() == sp)
        {
            retval.push_back(*i);
        }
    }

    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::list_particles_within_radius(
    const Real3& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dist(distance((*i).second.position(), pos) - (*i).second.radius());
        if (dist <= radius)
        {
            retval.push_back(std::make_pair(*i, dist));
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::list_particles_within_radius(
    const Real3& pos, const Real& radius, const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dist(distance((*i).second.position(), pos) - (*i).second.radius());
        if (dist <= radius)
        {
            if ((*i).first != ignore)
            {
                retval.push_back(std::make_pair(*i, dist));
            }
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::list_particles_within_radius(
    const Real3& pos, const Real& radius,
    const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dist(distance((*i).second.position(), pos) - (*i).second.radius());
        if (dist <= radius)
        {
            if ((*i).first != ignore1 && (*i).first != ignore2)
            {
                retval.push_back(std::make_pair(*i, dist));
            }
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

void ParticleSpaceVectorImpl::reset(const Real3& edge_lengths)
{
    base_type::t_ = 0.0;
    particles_.clear();
    index_map_.clear();

    for (Real3::size_type dim(0); dim < 3; ++dim)
    {
        if (edge_lengths[dim] <= 0)
        {
            throw std::invalid_argument("the edge length must be positive.");
        }
    }

    edge_lengths_ = edge_lengths;
}

} // ecell4
