#include "BDContainer2D.hpp"
#include <ecell4/core/Context.hpp>

namespace ecell4
{

namespace bd
{

Integer ParticleContainer2D::num_particles() const
{
    return particles_.size();
}

Integer ParticleContainer2D::num_particles(const Species& sp) const
{
    Integer retval = 0;
    SpeciesExpressionMatcher sexp(sp);
    for(species_to_particle_id_set_map_type::const_iterator
            iter(particle_pool_.begin()); iter != particle_pool_.end(); ++iter)
    {
        const Species target((*iter).first);
        if(sexp.match(target))
        {
            retval += (*iter).second.size();
        }
    }
    return retval;
}

Integer ParticleContainer2D::num_particles_exact(const Species& sp) const
{
    species_to_particle_id_set_map_type::const_iterator
        i(particle_pool_.find(sp.serial()));
    if (i == particle_pool_.end())
    {
        return 0;
    }
    return (*i).second.size();
}

Integer ParticleContainer2D::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for (species_to_particle_id_set_map_type::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
    {
        const Species tgt((*i).first);
        retval += sexp.count(tgt) * (*i).second.size();
    }
    return retval;
}

Integer ParticleContainer2D::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

std::vector<std::pair<ParticleID, Particle> >
ParticleContainer2D::list_particles() const
{
    return particles_;
}

std::vector<std::pair<ParticleID, Particle> >
ParticleContainer2D::list_particles(const Species& sp) const
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
ParticleContainer2D::list_particles_exact(const Species& sp) const
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


bool ParticleContainer2D::has_particle(const ParticleID& pid) const
{
    return (this->pid_to_pidx_.find(pid) != this->pid_to_pidx_.end());
}

bool
ParticleContainer2D::update_particle(
        const ParticleID& pid, const Particle& p, const face_id_type& fid)
{
    particle_container_type::iterator iter = this->find(pid);
    if(iter != particles_.end())
    {
        if(iter->second.species() != p.species())
        {
            particle_pool_[iter->second.species_serial()].erase(iter->first);
            particle_pool_[p.species_serial()].insert(pid);
        }
        this->update(std::make_pair(pid, p), fid);
        return false;
    }
    particle_pool_[p.species_serial()].insert(pid);
    this->update(std::make_pair(pid, p), fid);
    return true;
}

std::pair<ParticleID, Particle>
ParticleContainer2D::get_particle(const ParticleID& pid) const
{
    const particle_container_type::const_iterator iter = this->find(pid);
    if(iter == particles_.end()) throw NotFound("No such particle");
    return *iter;
}

void ParticleContainer2D::remove_particle(const ParticleID& pid)
{
    const std::pair<ParticleID, Particle> p = this->get_particle(pid);
    particle_pool_[p.second.species_serial()].erase(pid);
    this->erase(pid);
    return;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleContainer2D::list_particles_within_radius(
        const Real3& pos, const Real& radius) const
{
    const Real rad2 = radius * radius;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    for(particle_container_type::const_iterator
            iter = particles_.begin(); iter != particles_.end(); ++iter)
    {
        const Real3 dpos = pos - iter->second.position();
        const Real l2 = length_sq(dpos);

        const face_type& face = this->belonging_face(iter->first);
        // XXX: slow...
        const Real phi = std::acos(dot_product(dpos, face.normal()) / std::sqrt(l2));
        const Real theta = (phi < M_PI * 0.5) ? M_PI * 0.5 - phi : phi - M_PI * 0.5;
        const Real dist2 = l2 + rad2 - 2 * std::sqrt(l2 * rad2) * std::cos(theta);

        if(dist2 <= rad2)
            retval.push_back(std::make_pair(*iter, std::sqrt(dist2)));
    }
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleContainer2D::list_particles_within_radius(
        const Real3& pos, const Real& radius, const ParticleID& ignore) const
{
    const Real rad2 = radius * radius;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    for(particle_container_type::const_iterator
            iter = particles_.begin(); iter != particles_.end(); ++iter)
    {
        if(iter->first == ignore) continue;

        const Real3 dpos = pos - iter->second.position();
        const Real l2 = length_sq(dpos);

        const face_type& face = this->belonging_face(iter->first);
        // XXX: slow...
        const Real phi = std::acos(dot_product(dpos, face.normal()) / std::sqrt(l2));
        const Real theta = (phi < M_PI * 0.5) ? M_PI * 0.5 - phi : phi - M_PI * 0.5;
        const Real dist2 = l2 + rad2 - 2 * std::sqrt(l2 * rad2) * std::cos(theta);

        if(dist2 <= rad2)
            retval.push_back(std::make_pair(*iter, std::sqrt(dist2)));
    }
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleContainer2D::list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    const Real rad2 = radius * radius;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    for(particle_container_type::const_iterator
            iter = particles_.begin(); iter != particles_.end(); ++iter)
    {
        if(iter->first == ignore1 || iter->first == ignore2) continue;

        const Real3 dpos = pos - iter->second.position();
        const Real l2 = length_sq(dpos);

        const face_type& face = this->belonging_face(iter->first);
        // XXX: slow...
        const Real phi = std::acos(dot_product(dpos, face.normal()) / std::sqrt(l2));
        const Real theta = (phi < M_PI * 0.5) ? M_PI * 0.5 - phi : phi - M_PI * 0.5;
        const Real dist2 = l2 + rad2 - 2 * std::sqrt(l2 * rad2) * std::cos(theta);

        if(dist2 <= rad2)
            retval.push_back(std::make_pair(*iter, std::sqrt(dist2)));
    }
    return retval;
}


std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleContainer2D::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    const face_id_list& neighbors = polygon_.neighbor_faces(pos.second);
    for(face_id_list::const_iterator
            iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const particle_id_set& ps_on_face = this->particles_on_face(*iter);
        for(particle_id_set::const_iterator
                jter = ps_on_face.begin(); jter != ps_on_face.end(); ++jter)
        {
            const ParticleID pid = *jter;
            const face_id_type fid = this->belonging_faceid(pid);
            const particle_container_type::const_iterator piter =
                this->find(pid);
            const Real l = polygon_.distance(pos, std::make_pair(
                        piter->second.position(), fid)) - piter->second.radius();
            if(l <= radius)
                retval.push_back(std::make_pair(*piter, l));
        }
    }
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleContainer2D::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos,
        const Real& radius, const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    const face_id_list& neighbors = polygon_.neighbor_faces(pos.second);
    for(face_id_list::const_iterator
            iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const particle_id_set& ps_on_face = this->particles_on_face(*iter);
        for(particle_id_set::const_iterator
                jter = ps_on_face.begin(); jter != ps_on_face.end(); ++jter)
        {
            const ParticleID pid = *jter;
            if(pid == ignore) continue;
            const face_id_type fid = this->belonging_faceid(pid);
            const particle_container_type::const_iterator piter =
                this->find(pid);
            const Real l = polygon_.distance(pos, std::make_pair(
                        piter->second.position(), fid)) - piter->second.radius();
            if(l <= radius)
                retval.push_back(std::make_pair(*piter, l));
        }
    }
    return retval;    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleContainer2D::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    const face_id_list& neighbors = polygon_.neighbor_faces(pos.second);
    for(face_id_list::const_iterator
            iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const particle_id_set& ps_on_face = this->particles_on_face(*iter);
        for(particle_id_set::const_iterator
                jter = ps_on_face.begin(); jter != ps_on_face.end(); ++jter)
        {
            const ParticleID pid = *jter;
            if(pid == ignore1 || pid == ignore2) continue;
            const face_id_type fid = this->belonging_faceid(pid);
            const particle_container_type::const_iterator piter =
                this->find(pid);
            const Real l = polygon_.distance(pos, std::make_pair(
                        piter->second.position(), fid)) - piter->second.radius();
            if(l <= radius)
                retval.push_back(std::make_pair(*piter, l));
        }
    }    return retval;
}

std::pair<Real3, ParticleContainer2D::face_id_type>
ParticleContainer2D::apply_surface(
        const std::pair<Real3, face_id_type>& position,
        const Real3& displacement) const
{
    std::pair<std::pair<Real3, face_id_type>, Real3>
        state = std::make_pair(position, displacement);
    const Real cutoff = length_sq(displacement) * 1e-4;//tolerance 1e-8

    Integer retry_count = 100;
    while(--retry_count > 0 && length_sq(state.second) > cutoff)
    {
        state = polygon_.move_next_face(state.first, state.second);
    }
    if(retry_count == 0)
        std::cerr << "Warning: max move_next_face count exceeded" << std::endl;
//     else if(length_sq(state.second) != 0.0)
//         std::cerr << "Warning: displacement length cut." << std::endl;

    const face_id_type fid = state.first.second;
    Real3 newpos = state.first.first;

    const Triangle& face = this->polygon_.at(fid);
    const Real      dist = dot_product(face.normal(), newpos - face.vertex_at(0));
    if(std::abs(dist) > 1e-8)
        std::cerr << "Warning: particle floats over tolerance: " << dist << std::endl;
    newpos = newpos - face.normal() * dist;

    Barycentric<Real> bary(to_barycentric(newpos, face));
    if(!is_inside(bary, 1e-8))
        std::cerr << "Warning: particle is not in triangle: barycentric("
                  << bary << ")" << std::endl;
    bary = force_put_inside(bary);
    newpos = to_absolute(bary, face);

    return std::make_pair(newpos, fid);
}

void ParticleContainer2D::setup_polygon()
{
    this->polygon_.detect_connectivity();
    std::vector<face_id_type> face_ids(polygon_.list_faceids());
    for(std::vector<face_id_type>::const_iterator
            iter = face_ids.begin(); iter != face_ids.end(); ++iter)
    {
        particle_on_face_[*iter] = particle_id_set();
    }

    return;
}

}// bd
}// ecell4
