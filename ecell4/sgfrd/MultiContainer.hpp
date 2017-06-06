#ifndef ECELL4_SGFRD_MULTI_CONTAINER
#define ECELL4_SGFRD_MULTI_CONTAINER
#include <ecell4/sgfrd/SGFRDWorld.hpp>

namespace ecell4
{
namespace sgfrd
{

class MultiContainer
{
  public:
    typedef SGFRDWorld world_type;
    typedef world_type::polygon_type       polygon_type;
    typedef world_type::triangle_type      triangle_type;
    typedef world_type::face_id_type       face_id_type;
    typedef world_type::edge_id_type       edge_id_type;
    typedef world_type::vertex_id_type     vertex_id_type;
    typedef world_type::face_descripter    face_descripter;
    typedef world_type::edge_descripter    edge_descripter;
    typedef world_type::vertex_descripter  vertex_descripter;
    typedef world_type::local_index_type   local_index_type;
    typedef world_type::barycentric_type   barycentric_type;
    typedef world_type::molecule_info_type molecule_info_type;
    typedef world_type::structure_registrator_type structure_registrator_type;
    typedef world_type::particle_space_type        particle_space_type;
    typedef world_type::particle_container_type    particle_container_type;

  public:

    MultiContainer(world_type& w) : world_(w), registrator_(*(w.polygon())){}
    ~MultiContainer(){}

    void make_entry(const ParticleID& pid)
    {
        pcon_.push_back(world_.get_particle(pid));
        if(world_.is_on_face(pid))
            registrator_.emplace(pid, world_.get_face_id(pid));
        return;
    }

    particle_container_type&       list_particles()       {return pcon_;}
    particle_container_type const& list_particles() const {return pcon_;}

    std::size_t num_particles() const {return pcon_.size();}

    Real t() const {return world_.t();}

    face_id_type get_face_id(const ParticleID& pid) const
    {return world_.get_face_id(pid);}

    molecule_info_type get_molecule_info(const Species& sp) const
    {return world_.get_molecule_info(sp);}

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        return world_.update_particle(pid, p);
    }
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type fid)
    {
        if(registrator_.have(pid)) registrator_.update(pid, fid);
        else                       registrator_.emplace(pid, fid);
        return world_.update_particle(pid, p, fid);
    }

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p)
    {
        return world_.new_particle(p);
    }
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p, const face_id_type& fid)
    {
        const std::pair<std::pair<ParticleID, Particle>, bool> result =
            world_.new_particle(p, fid);
        if(result.second)
            registrator_.emplace(result.first.first, fid);
        return result;
    }

    void remove_particle(const ParticleID& pid)
    {
        if(registrator_.have(pid)) registrator_.remove(pid);
        return world_.remove_particle(pid);
    }
    void remove_particle(const ParticleID& pid, const face_id_type& fid)
    {
        registrator_.remove(pid, fid);
        return world_.remove_particle(pid);
    }

    // check particles associated with this Multi domain only.
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
            const std::pair<Real3, face_id_type>& pos, const Real& radius) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        const Real rad2 = radius * radius;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d2 = world_.distance_sq(std::make_pair(p.position(), fid), pos);
            if(d2 <= rad2) retval.push_back(std::make_pair(
                        std::make_pair(pid, p), std::sqrt(d2)));
        }
        std::sort(retval.begin(), retval.end(),
            ecell4::utils::pair_second_element_comparator<
                std::pair<ParticleID, Particle>, Real>());
        return retval;
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
            const std::pair<Real3, face_id_type>& pos, const Real& radius,
            const ParticleID& ignore) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        const Real rad2 = radius * radius;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore) continue;
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d2 = world_.distance_sq(std::make_pair(p.position(), fid), pos);
            if(d2 <= rad2) retval.push_back(std::make_pair(
                        std::make_pair(pid, p), std::sqrt(d2)));
        }
        std::sort(retval.begin(), retval.end(),
            ecell4::utils::pair_second_element_comparator<
                std::pair<ParticleID, Particle>, Real>());
        return retval;
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
            const std::pair<Real3, face_id_type>& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        const Real rad2 = radius * radius;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore1 || pid == ignore2) continue;
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d2 = world_.distance_sq(std::make_pair(p.position(), fid), pos);
            if(d2 <= rad2) retval.push_back(std::make_pair(
                        std::make_pair(pid, p), std::sqrt(d2)));
        }
        std::sort(retval.begin(), retval.end(),
            ecell4::utils::pair_second_element_comparator<
                std::pair<ParticleID, Particle>, Real>());
        return retval;
    }

    // return false if overlap exists.
    bool check_no_overlap(const std::pair<Real3, face_id_type>& pos,
            const Real& radius) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        const Real rad2 = radius * radius;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d2 = world_.distance_sq(std::make_pair(p.position(), fid), pos);
            if(d2 <= rad2) return false;
        }
        return true; // no overlap!
    }
    bool check_no_overlap(const std::pair<Real3, face_id_type>& pos,
            const Real& radius, const ParticleID& ignore) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        const Real rad2 = radius * radius;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore) continue;
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d2 = world_.distance_sq(std::make_pair(p.position(), fid), pos);
            if(d2 <= rad2) return false;
        }
        return true; // no overlap!
    }
    bool check_no_overlap(const std::pair<Real3, face_id_type>& pos,
            const Real& radius, const ParticleID& ignore1,
            const ParticleID& ignore2) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        const Real rad2 = radius * radius;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore1 || pid == ignore2) continue;
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d2 = world_.distance_sq(std::make_pair(p.position(), fid), pos);
            if(d2 <= rad2) return false;
        }
        return true; // no overlap!
    }


  private:

    world_type&                world_;
    structure_registrator_type registrator_;
    particle_container_type    pcon_;
};

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_MULTI_CONTAINER
