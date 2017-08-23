#include <ecell4/core/ParticleSpaceRTreeImpl.hpp>
#include <boost/geometry/algorithms/within.hpp>

namespace ecell4
{

void ParticleSpaceRTreeImpl::reset(const Real3& edge_lengths)
{
    base_type::t_ = 0.0;
    particles_.clear();
    idx_map_.clear();
    particle_pool_.clear();
    rtree_.clear();

    edge_lengths_ = edge_lengths_;
}

std::vector<Species> ParticleSpaceRTreeImpl::list_species() const
{
    std::vector<Species> retval;
    for(per_species_particle_id_set::const_iterator
        i(particle_pool_.begin()), e(particle_pool_.end()); i != e; ++i)
    {
        retval.push_back(Species(i->first));
    }
    return retval;
}

Integer ParticleSpaceRTreeImpl::num_particles(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for(per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()), e(particle_pool_.end()); i != e; ++i)
        if(sexp.match(i->first))
            retval += i->second.size();
    return retval;
}

Integer ParticleSpaceRTreeImpl::num_particles_exact(const Species& sp) const
{
    const per_species_particle_id_set::const_iterator
        i = particle_pool_.find(sp.serial());

    return (i == particle_pool_.end()) ? 0 : i->second.size();
}

Integer ParticleSpaceRTreeImpl::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for(per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()), e(particle_pool_.end()); i != e; ++i)
        retval += sexp.count(i->first) * i->second.size();
    return retval;
}

Integer ParticleSpaceRTreeImpl::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceRTreeImpl::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    SpeciesExpressionMatcher sexp(sp);

    for(particle_container_type::const_iterator
            i(particles_.begin()), e(particles_.end()); i != e; ++i)
        if(sexp.match(i->second.species()))
            retval.push_back(*i);

    return retval;
}
std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceRTreeImpl::list_particles_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    for(particle_container_type::const_iterator
            i(particles_.begin()), e(particles_.end()); i != e; ++i)
        if (i->second.species() == sp)
            retval.push_back(*i);

    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceRTreeImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    if(this->particles_.empty()) return retval;

    boost::container::static_vector<box_type, 8>
        boxes(1, self_type::make_box(pos, radius));

    const box_type boundary(Real3(0,0,0), this->edge_length_);
    if(!boost::geometry::within(boxes.front(), boundary))
    {// if the query box is out of periodic-boundary, split the query box
        split_box_by_boundary<0>(boxes);
        split_box_by_boundary<1>(boxes);
        split_box_by_boundary<2>(boxes);
    }

    for(boost::container::static_vector<box_type, 8>::const_iterator
            bxi(boxes.begin()), bxe(boxes.end()); bxi != bxe; ++bxi)
    {
        boost::container::small_vector<rtree_value_type, qsz> tmp;
        this->rtree_.query(boost::geometry::index::intersect(*bxi),
                           std::back_inserter(tmp));

        for(boost::container::small_vector<rtree_value_type, qsz>::const_iterator
                i(tmp.begin()), e(tmp.end()); i != e; ++i)
        {
            const ParticleID& pid = boost::get<1>(*i);
            const Particle&   p   = boost::get<2>(*i);
            const Real dist = length(p.position() - pos) - p.radius();

            if(dist < radius)
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), dist));
            }
        }
    }
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
list_particles_within_radius(const Real3& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    if(this->particles_.empty()) return retval;

    boost::container::static_vector<box_type, 8>
        boxes(1, self_type::make_box(pos, radius));

    const box_type boundary(Real3(0,0,0), this->edge_length_);
    if(!boost::geometry::within(boxes.front(), boundary))
    {// if the query box is out of periodic-boundary, split the query box
        split_box_by_boundary<0>(boxes);
        split_box_by_boundary<1>(boxes);
        split_box_by_boundary<2>(boxes);
    }

    for(boost::container::static_vector<box_type, 8>::const_iterator
            bxi(boxes.begin()), bxe(boxes.end()); bxi != bxe; ++bxi)
    {
        boost::container::small_vector<rtree_value_type, qsz> tmp;
        this->rtree_.query(boost::geometry::index::intersect(*bxi) &&
            boost::geometry::index::satisfies(particle_id_excluder(ignore)),
            std::back_inserter(tmp));

        for(boost::container::small_vector<rtree_value_type, qsz>::const_iterator
                i(tmp.begin()), e(tmp.end()); i != e; ++i)
        {
            const ParticleID& pid = boost::get<1>(*i);
            const Particle&   p   = boost::get<2>(*i);
            const Real dist = length(p.position() - pos) - p.radius();

            if(dist < radius)
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), dist));
            }
        }
    }
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
list_particles_within_radius(const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    if(this->particles_.empty()) return retval;

    boost::container::static_vector<box_type, 8>
        boxes(1, self_type::make_box(pos, radius));

    const box_type boundary(Real3(0,0,0), this->edge_length_);
    if(!boost::geometry::within(boxes.front(), boundary))
    {// if the query box is out of periodic-boundary, split the query box
        split_box_by_boundary<0>(boxes);
        split_box_by_boundary<1>(boxes);
        split_box_by_boundary<2>(boxes);
    }

    for(boost::container::static_vector<box_type, 8>::const_iterator
            bxi(boxes.begin()), bxe(boxes.end()); bxi != bxe; ++bxi)
    {
        boost::container::small_vector<rtree_value_type, qsz> tmp;
        this->rtree_.query(boost::geometry::index::intersect(*bxi) &&
            boost::geometry::index::satisfies(particle_id2_excluder(ignore1, ignore2)),
            std::back_inserter(tmp));

        for(boost::container::small_vector<rtree_value_type, qsz>::const_iterator
                i(tmp.begin()), e(tmp.end()); i != e; ++i)
        {
            const ParticleID& pid = boost::get<1>(*i);
            const Particle&   p   = boost::get<2>(*i);
            const Real dist = length(p.position() - pos) - p.radius();

            if(dist < radius)
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), dist));
            }
        }
    }
    return retval;
}

} // ecell4
