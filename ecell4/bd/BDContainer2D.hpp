#ifndef ECELL4_BD_CONTAINER_2D
#define ECELL4_BD_CONTAINER_2D
#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/exceptions.hpp>
#include "BDPolygon.hpp"
#include "map_const_at.hpp"
#include <set>

namespace ecell4
{

namespace bd
{

class ParticleContainer2D : public ParticleSpace
{
public:
    typedef ParticleSpace base_type;
    typedef base_type::particle_container_type particle_container_type;
    typedef particle_container_type::size_type container_index_type;
    typedef utils::get_mapper_mf<ParticleID, container_index_type>::type
            pid_to_particle_index_map_type;
    typedef std::set<ParticleID> particle_id_set;
    typedef std::map<Species::serial_type, particle_id_set>
            species_to_particle_id_set_map_type;

    typedef BDPolygon polygon_type;
    typedef polygon_type::triangle_type face_type;
    typedef polygon_type::face_id_type face_id_type;
    typedef std::vector<face_id_type> face_id_list;
    typedef utils::get_mapper_mf<ParticleID, face_id_type>::type
            pid_to_faceid_map_type;
    typedef utils::get_mapper_mf<face_id_type, particle_id_set>::type
            face_id_to_particle_id_set_map_type;

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
    list_particles_within_radius(const Real3& pos, const Real& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius,
        const ParticleID& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const;

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

    inline Real distance_sq(const Real3& pos1, const face_id_type& fid1,
                            const Real3& pos2, const face_id_type& fid2) const
    {
        return polygon_.distance_sq(
                std::make_pair(pos1, fid1), std::make_pair(pos2, fid2));
    }

    inline Real distance(const Real3& pos1, const face_id_type& fid1,
                         const Real3& pos2, const face_id_type& fid2) const
    {
        return polygon_.distance(
                std::make_pair(pos1, fid1), std::make_pair(pos2, fid2));
    }

    Real3 get_inter_position_vector(
            const std::pair<Real3, face_id_type>& lhs,
            const std::pair<Real3, face_id_type>& rhs) const
    {
        return polygon_.developed_direction(lhs, rhs);
    }

    bool has_particle(const ParticleID& pid) const;
    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        throw NotImplemented("2D::update_particle(pid, p)");
    }
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type& fid);

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const;
    void remove_particle(const ParticleID& pid);
    const particle_container_type& particles() const {return particles_;}

    // polygon
    std::pair<Real3, face_id_type>
    apply_surface(const std::pair<Real3, face_id_type>& position,
                  const Real3& displacement) const;

    polygon_type&       polygon()       {return polygon_;}
    polygon_type const& polygon() const {return polygon_;}

    void setup_polygon();

    face_type const& belonging_face(const ParticleID& pid) const
    {
        return polygon_.triangle_at(const_at(pid_to_fid_, pid));
    }

    face_id_type const& belonging_faceid(const ParticleID& pid) const
    {
        return const_at(pid_to_fid_, pid);
    }

    particle_id_set const& particles_on_face(const face_id_type& fid) const
    {
        return const_at(particle_on_face_, fid);
    }

#ifdef WITH_HDF5
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
#else
    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }
    virtual void load(const std::string& filename) const
    {
        throw NotSupported(
            "load(const std::string) is not supported by this space class");
    }
#endif

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

private:

    Real3                               edge_lengths_;
    polygon_type                        polygon_;
    pid_to_particle_index_map_type      pid_to_pidx_;
    pid_to_faceid_map_type              pid_to_fid_;
    species_to_particle_id_set_map_type particle_pool_;
    face_id_to_particle_id_set_map_type particle_on_face_;
    particle_container_type             particles_;
};

inline ParticleContainer2D::particle_container_type::iterator
ParticleContainer2D::find(const ParticleID& pid)
{
    const pid_to_particle_index_map_type::const_iterator
        iter(pid_to_pidx_.find(pid));
    if(iter == pid_to_pidx_.end()) return particles_.end();
    return particles_.begin() + iter->second;
}

inline ParticleContainer2D::particle_container_type::const_iterator
ParticleContainer2D::find(const ParticleID& pid) const
{
    const pid_to_particle_index_map_type::const_iterator
        iter(pid_to_pidx_.find(pid));
    if(iter == pid_to_pidx_.end()) return particles_.end();
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
        pid_to_pidx_[p.first] = idx;
        pid_to_fid_[p.first] = fid;
        particle_on_face_[fid].insert(p.first);
        return particles_.begin() + idx;
    }

    *old = p;
    const face_id_type old_face(pid_to_fid_[p.first]);

    if(old_face != fid)
    {
        particle_on_face_[old_face].erase(p.first);
        particle_on_face_[fid].insert(p.first);
        pid_to_fid_[p.first] = fid;
    }
    return old;
}

inline std::pair<ParticleContainer2D::particle_container_type::iterator, bool>
ParticleContainer2D::update(
        const std::pair<ParticleID, Particle>& pid, const face_id_type fid)
{
    const pid_to_particle_index_map_type::const_iterator
        iter(this->pid_to_pidx_.find(pid.first));
    if(iter != pid_to_pidx_.end())
        return std::make_pair(
                this->update(particles_.begin() + iter->second, pid, fid),
                false);

    return std::make_pair(this->update(particles_.end(), pid, fid), true);
}

inline bool
ParticleContainer2D::erase(const particle_container_type::iterator& iter)
{
    if(particles_.end() == iter)
        return false;

    const container_index_type old_idx(std::distance(particles_.begin(), iter));
    const face_id_type old_face(this->pid_to_fid_[iter->first]);

    this->pid_to_pidx_.erase(iter->first);
    this->pid_to_fid_.erase(iter->first);
    this->particle_on_face_[old_face].erase(iter->first);

    const container_index_type last_idx(particles_.size() - 1);
    if(old_idx < last_idx)
    {
        // exchange the last element and the element to remove
        // because pop_back is efficient
        const std::pair<ParticleID, Particle>& last(particles_[last_idx]);
        pid_to_pidx_[last.first] = old_idx;
        *iter = last;
    }
    particles_.pop_back();

    return true;
}

inline bool ParticleContainer2D::erase(const ParticleID& pid)
{
    const pid_to_particle_index_map_type::const_iterator
        iter(this->pid_to_pidx_.find(pid));
    if(this->pid_to_pidx_.end() == iter) return false;
    return this->erase(particles_.begin() + iter->second);
}

}// bd
}// ecell4

#endif /* ECELL4_BD_BD_CONTAINER_2D */
