#include <ecell4/core/ParticleSpaceRTreeImpl.hpp>
#include <ecell4/core/comparators.hpp>
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

    this->edge_lengths_ = edge_lengths;
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
    {
        const Species target(i->first);
        if(sexp.match(target))
        {
            retval += i->second.size();
        }
    }
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
    {
        const Species target(i->first);
        retval += sexp.count(target) * i->second.size();
    }
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
    {
        if(sexp.match(i->second.species()))
        {
            retval.push_back(*i);
        }
    }

    return retval;
}
std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceRTreeImpl::list_particles_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    for(particle_container_type::const_iterator
            i(particles_.begin()), e(particles_.end()); i != e; ++i)
    {
        if (i->second.species() == sp)
        {
            retval.push_back(*i);
        }
    }

    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceRTreeImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    if(this->particles_.empty())
    {
        return retval;
    }

    query_boxes_container_type
        boxes(1, self_type::make_box(pos, radius+this->max_radius_));

    const box_type boundary(Real3(0,0,0), this->edge_lengths_);
    if(not boost::geometry::within(boxes.front(), boundary))
    {// if the query box is out of periodic-boundary, split the query box
        this->split_box_by_boundary<0>(boxes);
        this->split_box_by_boundary<1>(boxes);
        this->split_box_by_boundary<2>(boxes);
    }

    for(query_boxes_container_type::const_iterator
            bxi(boxes.begin()), bxe(boxes.end()); bxi != bxe; ++bxi)
    {
        query_result_container_type tmp;
        this->rtree_.query(boost::geometry::index::intersects(*bxi),
                           std::back_inserter(tmp));

        for(query_result_container_type::const_iterator
                i(tmp.begin()), e(tmp.end()); i != e; ++i)
        {
            const ParticleID& pid = boost::get<1>(*i);
            const Particle&   p   = boost::get<2>(*i);
            const Real dist = this->distance(p.position(), pos) - p.radius();

            if(dist <= radius)
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(), utils::pair_second_element_comparator<
            std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceRTreeImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius, const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    if(this->particles_.empty())
    {
        return retval;
    }

    query_boxes_container_type
        boxes(1, self_type::make_box(pos, radius+this->max_radius_));

    const box_type boundary(Real3(0,0,0), this->edge_lengths_);
    if(not boost::geometry::within(boxes.front(), boundary))
    {// if the query box is out of periodic-boundary, split the query box
        this->split_box_by_boundary<0>(boxes);
        this->split_box_by_boundary<1>(boxes);
        this->split_box_by_boundary<2>(boxes);
    }

    for(query_boxes_container_type::const_iterator
            bxi(boxes.begin()), bxe(boxes.end()); bxi != bxe; ++bxi)
    {
        query_result_container_type tmp;
        this->rtree_.query(boost::geometry::index::intersects(*bxi) &&
            boost::geometry::index::satisfies(particle_id_excluder(ignore)),
            std::back_inserter(tmp));

        for(query_result_container_type::const_iterator
                i(tmp.begin()), e(tmp.end()); i != e; ++i)
        {
            const ParticleID& pid = boost::get<1>(*i);
            const Particle&   p   = boost::get<2>(*i);
            const Real dist = this->distance(p.position(), pos) - p.radius();

            BOOST_ASSERT(pid != ignore);

            if(dist <= radius)
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(), utils::pair_second_element_comparator<
            std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceRTreeImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    if(this->particles_.empty())
    {
        return retval;
    }

    query_boxes_container_type
        boxes(1, self_type::make_box(pos, radius+this->max_radius_));

    const box_type boundary(Real3(0,0,0), this->edge_lengths_);

    if(not boost::geometry::within(boxes.front(), boundary))
    {// if the query box is out of periodic-boundary, split the query box
        this->split_box_by_boundary<0>(boxes);
        this->split_box_by_boundary<1>(boxes);
        this->split_box_by_boundary<2>(boxes);
    }

    for(query_boxes_container_type::const_iterator
            bxi(boxes.begin()), bxe(boxes.end()); bxi != bxe; ++bxi)
    {
        query_result_container_type tmp;
        this->rtree_.query(boost::geometry::index::intersects(*bxi) &&
            boost::geometry::index::satisfies(particle_id2_excluder(ignore1, ignore2)),
            std::back_inserter(tmp));

        for(query_result_container_type::const_iterator
                i(tmp.begin()), e(tmp.end()); i != e; ++i)
        {
            const ParticleID& pid = boost::get<1>(*i);
            const Particle&   p   = boost::get<2>(*i);
            const Real dist = this->distance(p.position(), pos) - p.radius();

            BOOST_ASSERT(pid != ignore1);
            BOOST_ASSERT(pid != ignore2);

            if(dist <= radius)
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(), utils::pair_second_element_comparator<
            std::pair<ParticleID, Particle>, Real>());

    return retval;
}

} // ecell4
