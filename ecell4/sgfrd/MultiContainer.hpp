#ifndef ECELL4_SGFRD_MULTI_CONTAINER
#define ECELL4_SGFRD_MULTI_CONTAINER
#include <ecell4/sgfrd/SGFRDWorld.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

namespace ecell4
{
namespace sgfrd
{

class MultiContainer
{
  public:
    typedef SGFRDWorld world_type;
    typedef Polygon  polygon_type;
    typedef world_type::FaceID   FaceID;
    typedef world_type::EdgeID   EdgeID;
    typedef world_type::VertexID VertexID;
    typedef world_type::structure_registrator_type structure_registrator_type;
    typedef world_type::particle_space_type        particle_space_type;
    typedef world_type::particle_container_type    particle_container_type;

  public:

    MultiContainer(world_type& w) : world_(w), registrator_(*(w.polygon())){}
    ~MultiContainer(){}

    bool make_entry(const ParticleID& pid)
    {
        if(this->find_(pid) != pcon_.end()) return false;

        pcon_.push_back(world_.get_particle(pid));
        if(world_.is_on_face(pid))
            registrator_.emplace(pid, world_.get_face_id(pid));
        return true;
    }

    particle_container_type&       list_particles()       {return pcon_;}
    particle_container_type const& list_particles() const {return pcon_;}

    std::size_t num_particles() const {return pcon_.size();}

    Real t() const {return world_.t();}

    FaceID get_face_id(const ParticleID& pid) const
    {return world_.get_face_id(pid);}

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        *(this->find_(pid)) = std::make_pair(pid, p);
        return world_.update_particle(pid, p);
    }
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const FaceID fid)
    {
        if(registrator_.have(pid))
        {
            registrator_.update(pid, fid);
        }
        else
        {
            registrator_.emplace(pid, fid);
        }
        *(this->find_(pid)) = std::make_pair(pid, p);

        return world_.update_particle(pid, p, fid);
    }

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p)
    {
        const std::pair<std::pair<ParticleID, Particle>, bool> result =
            world_.new_particle(p);
        if(result.second)
        {
            this->pcon_.push_back(result.first);
        }
        return result;
    }
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p, const FaceID& fid)
    {
        const std::pair<std::pair<ParticleID, Particle>, bool> result =
            world_.new_particle(p, fid);
        if(result.second)
        {
            this->pcon_.push_back(result.first);
            registrator_.emplace(result.first.first, fid);
        }
        return result;
    }

    void remove_particle(const ParticleID& pid)
    {
        if(registrator_.have(pid))
        {
            registrator_.remove(pid);
        }
        const particle_container_type::iterator to_be_removed = this->find_(pid);
        this->pcon_.erase(to_be_removed);

        return world_.remove_particle(pid);
    }
    void remove_particle(const ParticleID& pid, const FaceID& fid)
    {
        registrator_.remove(pid, fid);
        const particle_container_type::iterator to_be_removed = this->find_(pid);
        this->pcon_.erase(to_be_removed);
        return world_.remove_particle(pid);
    }

    // check particles associated with this Multi domain only.
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
            const std::pair<Real3, FaceID>& pos, const Real& radius) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d = world_.distance(std::make_pair(p.position(), fid), pos);
            if(d <= (radius + p.radius()))
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), d));
            }
        }
        std::sort(retval.begin(), retval.end(),
            ecell4::utils::pair_second_element_comparator<
                std::pair<ParticleID, Particle>, Real>());
        return retval;
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
            const std::pair<Real3, FaceID>& pos, const Real& radius,
            const ParticleID& ignore) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore){continue;}
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d = world_.distance(std::make_pair(p.position(), fid), pos);
            if(d <= (radius + p.radius()))
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), d));
            }
        }
        std::sort(retval.begin(), retval.end(),
            ecell4::utils::pair_second_element_comparator<
                std::pair<ParticleID, Particle>, Real>());
        return retval;
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
            const std::pair<Real3, FaceID>& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore1 || pid == ignore2){continue;}
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real d = world_.distance(std::make_pair(p.position(), fid), pos);
            if(d <= (radius + p.radius()))
            {
                retval.push_back(std::make_pair(std::make_pair(pid, p), d));
            }
        }
        std::sort(retval.begin(), retval.end(),
            ecell4::utils::pair_second_element_comparator<
                std::pair<ParticleID, Particle>, Real>());
        return retval;
    }

    // return false if overlap exists.
    bool check_no_overlap(const std::pair<Real3, FaceID>& pos,
            const Real& radius) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real dist =
                this->world_.distance(std::make_pair(p.position(), fid), pos);
            if(dist <= radius + p.radius())
            {
                return false; // overlaps!
            }
        }
        return true; // no overlap!
    }
    bool check_no_overlap(const std::pair<Real3, FaceID>& pos,
            const Real& radius, const ParticleID& ignore) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore){continue;}

            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real dist =
                this->world_.distance(std::make_pair(p.position(), fid), pos);
            if(dist <= radius + p.radius())
            {
                return false; // overlaps!
            }
        }
        return true; // no overlap!
    }
    bool check_no_overlap(const std::pair<Real3, FaceID>& pos,
            const Real& radius, const ParticleID& ignore1,
            const ParticleID& ignore2) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), pcon_)
        {
            if(pid == ignore1 || pid == ignore2){continue;}

            BOOST_AUTO(fid, registrator_.structure_on(pid));
            const Real dist =
                this->world_.distance(std::make_pair(p.position(), fid), pos);
            if(dist <= radius + p.radius())
            {
                return false; // overlaps!
            }
        }
        return true; // no overlap!
    }

    world_type const& world() const throw() {return world_;}

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {
        return this->world_.get_particle(pid);
    }

  private:

    particle_container_type::iterator find_(const ParticleID& pid)
    {
        return std::find_if(pcon_.begin(), pcon_.end(),
                  ecell4::utils::pair_first_element_unary_predicator<
                      ParticleID, Particle>(pid));
    }

  private:

    world_type&                world_;
    structure_registrator_type registrator_;
    particle_container_type    pcon_;
};

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_MULTI_CONTAINER
