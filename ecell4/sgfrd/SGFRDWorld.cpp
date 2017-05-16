#include <ecell4/sgfrd/SGFRDWorld.hpp>

namespace ecell4
{
namespace sgfrd
{

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::new_particle(const Particle& p)
{
    const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlap3 = list_particles_within_radius(p.position(), p.radius());
    if(!overlap3.empty())
    {
        return std::make_pair(std::make_pair(pidgen_(), p), false);
    }
    const ParticleID pid = pidgen_();
    return std::make_pair(std::make_pair(pid, p), update_particle(pid, p));
}

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::new_particle(const Particle& p, const face_id_type& fid)
{
    // XXX: consider particle shape and split overlap check for 2d and 3d
    const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlap3 = list_particles_within_radius(p.position(), p.radius());
    if(!overlap3.empty())
    {
        return std::make_pair(std::make_pair(pidgen_(), p), false);
    }
    const ParticleID pid = pidgen_();
    return std::make_pair(std::make_pair(pid, p), update_particle(pid, p, fid));
}


std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) retval.push_back(std::make_pair(pp, dist));
        }
    }

    std::vector<face_id_type> const& neighbors =
        polygon_->neighbor_faces(pos.second);
    for(std::vector<face_id_type>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = registrator_.elements_over(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) retval.push_back(std::make_pair(pp, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) continue;
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) retval.push_back(std::make_pair(pp, dist));
        }
    }

    std::vector<face_id_type> const& neighbors =
        polygon_->neighbor_faces(pos.second);
    for(std::vector<face_id_type>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) continue;
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) retval.push_back(std::make_pair(pp, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) continue;
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) retval.push_back(std::make_pair(pp, dist));
        }
    }

    std::vector<face_id_type> const& neighbors =
        polygon_->neighbor_faces(pos.second);
    for(std::vector<face_id_type>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) continue;
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) retval.push_back(std::make_pair(pp, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

}// sgfrd
}// ecell4
