#include "ParticleSpaceCellListImpl.hpp"
#include "Context.hpp"
#include "comparators.hpp"


namespace ecell4
{

void ParticleSpaceCellListImpl::reset(const Real3& edge_lengths)
{
    base_type::t_ = 0.0;
    particles_.clear();
    rmap_.clear();
    particle_pool_.clear();

    for (matrix_type::size_type i(0); i < matrix_.shape()[0]; ++i)
    {
        for (matrix_type::size_type j(0); j < matrix_.shape()[1]; ++j)
        {
            for (matrix_type::size_type k(0); k < matrix_.shape()[2]; ++k)
            {
                matrix_[i][j][k].clear();
            }
        }
    }

    for (Real3::size_type dim(0); dim < 3; ++dim)
    {
        if (edge_lengths[dim] <= 0)
        {
            throw std::invalid_argument("the edge length must be positive.");
        }
    }

    edge_lengths_ = edge_lengths;
    // throw NotImplemented("Not implemented yet.");
}

bool ParticleSpaceCellListImpl::update_particle(
    const ParticleID& pid, const Particle& p)
{
    particle_container_type::iterator i(find(pid));
    if (i != particles_.end())
    {
        if ((*i).second.species() != p.species())
        {
            particle_pool_[(*i).second.species_serial()].erase((*i).first);
            particle_pool_[p.species_serial()].insert(pid);
        }
        this->update(i, std::make_pair(pid, p));
        return false;
    }

    this->update(std::make_pair(pid, p));
    // const bool succeeded(this->update(std::make_pair(pid, p)).second);
    // BOOST_ASSERT(succeeded);

    particle_pool_[p.species_serial()].insert(pid);
    return true;
}

std::pair<ParticleID, Particle> ParticleSpaceCellListImpl::get_particle(
    const ParticleID& pid) const
{
    particle_container_type::const_iterator i(this->find(pid));
    if (i == particles_.end())
    {
        throw NotFound("No such particle.");
    }
    return (*i);
}

bool ParticleSpaceCellListImpl::has_particle(const ParticleID& pid) const
{
    return (this->find(pid) != particles_.end());
}

void ParticleSpaceCellListImpl::remove_particle(const ParticleID& pid)
{
    //XXX: In contrast to the original ParticleContainer in epdp,
    //XXX: this remove_particle throws an error when no corresponding
    //XXX: particle is found.
    std::pair<ParticleID, Particle> pp(get_particle(pid)); //XXX: may raise an error.
    particle_pool_[pp.second.species_serial()].erase(pid);
    this->erase(pid);
}

Integer ParticleSpaceCellListImpl::num_particles() const
{
    return particles_.size();
}

Integer ParticleSpaceCellListImpl::num_particles(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for (per_species_particle_id_set::const_iterator i(particle_pool_.begin());
        i != particle_pool_.end(); ++i)
    {
        const Species tgt((*i).first);
        if (sexp.match(tgt))
        {
            retval += (*i).second.size();
        }
    }
    return retval;
}

Integer ParticleSpaceCellListImpl::num_particles_exact(const Species& sp) const
{
    per_species_particle_id_set::const_iterator i(particle_pool_.find(sp.serial()));
    if (i == particle_pool_.end())
    {
        return 0;
    }
    return (*i).second.size();
}

Integer ParticleSpaceCellListImpl::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for (per_species_particle_id_set::const_iterator i(particle_pool_.begin());
        i != particle_pool_.end(); ++i)
    {
        const Species tgt((*i).first);
        retval += sexp.count(tgt) * (*i).second.size();
    }
    return retval;
}

Integer ParticleSpaceCellListImpl::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

std::vector<std::pair<ParticleID, Particle> >
    ParticleSpaceCellListImpl::list_particles() const
{
    return particles_;
}

std::vector<std::pair<ParticleID, Particle> >
    ParticleSpaceCellListImpl::list_particles(const Species& sp) const
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
    ParticleSpaceCellListImpl::list_particles_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    // per_species_particle_id_set::const_iterator
    //     i(particle_pool_.find(sp.serial()));
    // if (i == particle_pool_.end())
    // {
    //     //XXX: In the original, this raises an error,
    //     //XXX: but returns an empty vector here.
    //     return retval;
    // }
    // retval.reserve((*i).second.size());

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
    ParticleSpaceCellListImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    // MatrixSpace::each_neighbor_cyclic
    if (particles_.size() == 0)
    {
        return retval;
    }

    cell_index_type idx(this->index(pos));

    // MatrixSpace::each_neighbor_cyclic_loops
    cell_offset_type off;
    for (off[2] = -1; off[2] <= 1; ++off[2])
    {
        for (off[1] = -1; off[1] <= 1; ++off[1])
        {
            for (off[0] = -1; off[0] <= 1; ++off[0])
            {
                cell_index_type newidx(idx);
                const Real3 stride(this->offset_index_cyclic(newidx, off));
                const cell_type& c(this->cell(newidx));
                for (cell_type::const_iterator i(c.begin()); i != c.end(); ++i)
                {
                    // neighbor_filter::operator()
                    particle_container_type::const_iterator
                        itr(particles_.begin() + (*i));
                    // particle_container_type::const_iterator itr = particles_.begin();
                    // std::advance(itr, *i);

                    const Real dist(
                        length((*itr).second.position() + stride - pos)
                        - (*itr).second.radius());
                    if (dist < radius)
                    {
                        // overlap_checker::operator()
                        retval.push_back(
                            std::make_pair(*itr, dist));
                    }
                }
            }
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    ParticleSpaceCellListImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    // MatrixSpace::each_neighbor_cyclic
    if (particles_.size() == 0)
    {
        return retval;
    }

    cell_index_type idx(this->index(pos));

    // MatrixSpace::each_neighbor_cyclic_loops
    cell_offset_type off;
    for (off[2] = -1; off[2] <= 1; ++off[2])
    {
        for (off[1] = -1; off[1] <= 1; ++off[1])
        {
            for (off[0] = -1; off[0] <= 1; ++off[0])
            {
                cell_index_type newidx(idx);
                const Real3 stride(this->offset_index_cyclic(newidx, off));
                const cell_type& c(this->cell(newidx));
                for (cell_type::const_iterator i(c.begin()); i != c.end(); ++i)
                {
                    // neighbor_filter::operator()
                    particle_container_type::const_iterator
                        itr(particles_.begin() + (*i));

                    const Real dist(
                        length((*itr).second.position() + stride - pos)
                        - (*itr).second.radius());
                    if (dist < radius)
                    {
                        // overlap_checker::operator()
                        if ((*itr).first != ignore)
                        {
                            retval.push_back(
                                std::make_pair(*itr, dist));
                        }
                    }
                }
            }
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    ParticleSpaceCellListImpl::list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    // MatrixSpace::each_neighbor_cyclic
    if (particles_.size() == 0)
    {
        return retval;
    }

    cell_index_type idx(this->index(pos));

    // MatrixSpace::each_neighbor_cyclic_loops
    cell_offset_type off;
    for (off[2] = -1; off[2] <= 1; ++off[2])
    {
        for (off[1] = -1; off[1] <= 1; ++off[1])
        {
            for (off[0] = -1; off[0] <= 1; ++off[0])
            {
                cell_index_type newidx(idx);
                const Real3 stride(this->offset_index_cyclic(newidx, off));
                const cell_type& c(this->cell(newidx));
                for (cell_type::const_iterator i(c.begin()); i != c.end(); ++i)
                {
                    // neighbor_filter::operator()
                    particle_container_type::const_iterator
                        itr(particles_.begin() + (*i));

                    const Real dist(
                        length((*itr).second.position() + stride - pos)
                        - (*itr).second.radius());
                    if (dist < radius)
                    {
                        // overlap_checker::operator()
                        if ((*itr).first != ignore1 && (*itr).first != ignore2)
                        {
                            retval.push_back(
                                std::make_pair(*itr, dist));
                        }
                    }
                }
            }
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

};
