#ifndef ECELL4_BD_BD_CONTAINER_2D
#define ECELL4_BD_BD_CONTAINER_2D

#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/exceptions.hpp>
#include "BDPolygon.hpp"
#include <set>

namespace ecell4
{

namespace bd
{

class ParticleContainer2D : public ParticleSpace
{
public:
    typedef ParticleSpace base_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef typename particle_container_type::size_type container_index_type;
    typedef utils::get_mapper_mf<ParticleID, container_index_type>::type 
        pid_to_particle_index_type;
    typedef std::set<ParticleID> particle_id_set;
    typedef std::map<Species::serial_type, particle_id_set>
        per_species_particle_id_set;

    typedef BDPolygon polygon_type;
    typedef typename polygon_type::face_type face_type;
    typedef typename polygon_type::face_id_type face_id_type;
    typedef utils::get_mapper_mf<ParticleID, face_id_type>::type 
        pid_to_faceid_type;
    typedef utils::get_mapper_mf<face_id_type, particle_id_set>::type
        per_faces_particle_id_set;

public:

    ParticleContainer2D(){}
    ParticleContainer2D(const Real3& edge_lengths): edge_lengths_(edge_lengths){}
    ~ParticleContainer2D(){}

    const Real3& edge_lengths() const {return edge_lengths_;}
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    std::vector<std::pair<ParticleID, Particle> >
        list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles_exact(const Species& sp) const;
 
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius) const
    {
        throw NotImplemented("2D");
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius,
        const ParticleID& ignore) const
    {
        throw NotImplemented("2D");
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        throw NotImplemented("2D");
    }
 
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const std::pair<Real3, face_id_type>& pos, const Real& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const std::pair<Real3, face_id_type>& pos, const Real& radius,
            const ParticleID& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const std::pair<Real3, face_id_type>& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const;

    bool has_particle(const ParticleID& pid) const;
    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        throw NotImplemented("2D");
    }
    bool update_particle(const ParticleID& pid, const Particle& p, const face_id_type& fid);
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const;
    void remove_particle(const ParticleID& pid);
    const particle_container_type& particles() const {return particles_;}

    // polygon
    std::pair<Real3, face_id_type>
    apply_surface(const std::pair<Real3, face_id_type>& position, 
                  const Real3& displacement) const;


    polygon_type&       polygon()       {return polygon_;}
    polygon_type const& polygon() const {return polygon_;}

    face_type const& face_on(const ParticleID& pid)
    {
        return polygon_.at(fmap_[pid]);
    }

    face_id_type const& face_id_on(const ParticleID& pid)
    {
        return fmap_[pid];
    }

    void save_hdf5(H5::Group* root) const
    {
        throw NotSupported("2D hdf");
    }
    void load_hdf5(const H5::Group& root)
    {
        throw NotSupported("2D hdf");
    }
    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }

private:

    particle_container_type::iterator       find(const ParticleID& pid);
    particle_container_type::const_iterator find(const ParticleID& pid) const;

    particle_container_type::iterator
    update(const particle_container_type::iterator& old,
           const std::pair<ParticleID, Particle>& pid,
           const face_id_type fid);

    std::pair<particle_container_type::iterator, bool>
    update(const std::pair<ParticleID, Particle>& pid,
           const face_id_type fid);

    bool erase(const particle_container_type::iterator& pid);
    bool erase(const ParticleID& pid);

    struct particle_finder
        : public std::unary_function<std::pair<ParticleID, face_id_type>, bool>
    {
        particle_finder(const ParticleID& pid): pid_(pid){}

        bool operator()(std::pair<ParticleID, face_id_type> ptof) const
        {
            return ptof.first == pid_;
        }
      protected:
        ParticleID pid_;
    };


private:

    Real3                       edge_lengths_;
    polygon_type                polygon_;
    pid_to_particle_index_type  rmap_;
    pid_to_faceid_type          fmap_;
    per_species_particle_id_set particle_pool_;
    per_faces_particle_id_set   particle_face_;
    particle_container_type     particles_;
};

inline ParticleContainer2D::particle_container_type::iterator
ParticleContainer2D::find(const ParticleID& pid)
{
    const pid_to_particle_index_type::const_iterator iter = rmap_.find(pid);
    if(iter == rmap_.end()) return particles_.end();
    return particles_.begin() + iter->second;
}

inline ParticleContainer2D::particle_container_type::const_iterator
ParticleContainer2D::find(const ParticleID& pid) const
{
    const pid_to_particle_index_type::const_iterator iter = rmap_.find(pid);
    if(iter == rmap_.end()) return particles_.end();
    return particles_.begin() + iter->second;
}

inline ParticleContainer2D::particle_container_type::iterator
ParticleContainer2D::update(const particle_container_type::iterator& old,
                            const std::pair<ParticleID, Particle>& p,
                            const face_id_type fid)
{
    if(old == particles_.end())
    {
        const container_index_type idx = particles_.size();
        particles_.push_back(p);
        rmap_[p.first] = idx;
        fmap_[p.first] = fid;
        particle_face_[fid].insert(p.first);
        return particles_.begin() + idx;
    }

    *old = p;
    const face_id_type old_face = fmap_[p.first];

    if(old_face != fid)
    {
        particle_face_[old_face].erase(p.first);
        particle_face_[fid].insert(p.first);
        fmap_[p.first] = fid;
    }
    return old;
}

inline std::pair<ParticleContainer2D::particle_container_type::iterator, bool>
ParticleContainer2D::update(
        const std::pair<ParticleID, Particle>& pid, const face_id_type fid)
{
    pid_to_particle_index_type::const_iterator iter = this->rmap_.find(pid.first);
    if(iter != rmap_.end())
        return std::make_pair(this->update(particles_.begin() + iter->second, pid, fid), false);

    return std::make_pair(this->update(particles_.end(), pid, fid), true);
}

inline bool
ParticleContainer2D::erase(const particle_container_type::iterator& iter)
{
    if(particles_.end() == iter)
        return false;

    const container_index_type old_idx = std::distance(particles_.begin(), iter);
    const face_id_type old_face = this->fmap_[iter->first];

    this->rmap_.erase(iter->first);
    this->fmap_.erase(iter->first);
    this->particle_face_[old_face].erase(iter->first);

    const container_index_type last_idx = particles_.size() - 1;
    if(old_idx < last_idx)
    {// exchange the last element and the element to remove
        const std::pair<ParticleID, Particle>& last = particles_[last_idx];
        rmap_[last.first] = old_idx;
        *iter = last;
    }
    particles_.pop_back();

    return true;
}

inline bool ParticleContainer2D::erase(const ParticleID& pid)
{
    const pid_to_particle_index_type::const_iterator iter = this->rmap_.find(pid);
    if(this->rmap_.end() == iter) return false;
    return this->erase(particles_.begin() + iter->second);
}

}// bd
}// ecell4

#endif /* ECELL4_BD_BD_CONTAINER_2D */
