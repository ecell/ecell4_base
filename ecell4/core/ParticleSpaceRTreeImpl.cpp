#include <ecell4/core/ParticleSpaceRTreeImpl.hpp>
#include <ecell4/core/comparators.hpp>

namespace ecell4
{

void ParticleSpaceRTreeImpl::reset(const Real3& edge_lengths)
{
    this->t_ = 0.0;
    this->particle_pool_.clear();
    this->rtree_.clear();
    this->rtree_.edge_lengths() = edge_lengths;
    return;
}

std::vector<Species> ParticleSpaceRTreeImpl::list_species() const
{
    std::vector<Species> retval;
    for (const auto& pidp : rtree_.list_objects())
    {
        const Species& sp(pidp.second.species());
        if(std::find(retval.begin(), retval.end(), sp) == retval.end())
        {
            retval.push_back(sp);
        }
    }
    return retval;
}

Integer ParticleSpaceRTreeImpl::num_particles(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for(const auto& idset : particle_pool_)
    {
        const Species target(idset.first);
        if(sexp.match(target))
        {
            retval += idset.second.size();
        }
    }
    return retval;
}

Integer ParticleSpaceRTreeImpl::num_particles_exact(const Species& sp) const
{
    const auto i = particle_pool_.find(sp.serial());
    return (i == particle_pool_.end()) ? 0 : i->second.size();
}

Integer ParticleSpaceRTreeImpl::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for(const auto& idset : particle_pool_)
    {
        const Species target(idset.first);
        retval += sexp.count(target) * idset.second.size();
    }
    return retval;
}

Integer ParticleSpaceRTreeImpl::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

std::vector<std::pair<ParticleID, Particle>>
ParticleSpaceRTreeImpl::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle>> retval;
    SpeciesExpressionMatcher sexp(sp);

    for(const auto& pidp : rtree_.list_objects())
    {
        if(sexp.match(pidp.second.species()))
        {
            retval.push_back(pidp);
        }
    }
    return retval;
}
std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceRTreeImpl::list_particles_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle>> retval;

    for(const auto& pidp : rtree_.list_objects())
    {
        if (pidp.second.species() == sp)
        {
            retval.push_back(pidp);
        }
    }
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
ParticleSpaceRTreeImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> list;
    this->query_impl(make_intersection_query(pos, radius,
        [](const value_type&) noexcept -> bool {
            return false;
        }), std::back_inserter(list));

    return list;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
ParticleSpaceRTreeImpl::list_particles_within_radius(
    const Real3& pos, const Real& radius, const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> list;

    this->query_impl(make_intersection_query(pos, radius,
        [&ignore](const value_type& pidp) noexcept -> bool {
            return pidp.first == ignore;
        }), std::back_inserter(list));

    return list;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
ParticleSpaceRTreeImpl::list_particles_within_radius(
    const Real3& pos, const Real& radius, const ParticleID& ignore1,
    const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> list;

    this->query_impl(make_intersection_query(pos, radius,
        [&ignore1, &ignore2](const value_type& pidp) noexcept -> bool {
            return pidp.first == ignore1 || pidp.first == ignore2;
        }), std::back_inserter(list));

    return list;
}

} // ecell4
