#ifndef ECELL4_PARTICLE_SPACE_POLYGON_IMPL
#define ECELL4_PARTICLE_SPACE_POLYGON_IMPL
#include <ecell4/core/converter.hpp>
#include <functional>
#include <ecell4/core/comparators.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/sgfrd/StructuralContainer.hpp>
#include <set>

namespace ecell4
{
namespace sgfrd
{

template<typename T_traits>
class ParticleSpacePolygonImpl : public ParticleSpace
{
public:

    typedef ParticleSpace base_type;
    typedef ParticleSpace::particle_container_type particle_container_type;
    typedef typename ecell4::utils::get_mapper_mf<
        ParticleID, particle_container_type::size_type>::type pid_to_idx_map_type;

    typedef std::set<ParticleID> particle_id_set;
    typedef std::map<Species::serial_type, particle_id_set>
        per_species_particle_id_set;

    typedef T_traits traits_type;
    typedef Polygon<traits_type> polygon_type;
    typedef typename polygon_type::triangle_type     triangle_type;
    typedef typename polygon_type::face_id_type      face_id_type;
    typedef typename polygon_type::edge_id_type      edge_id_type;
    typedef typename polygon_type::vertex_id_type    vertex_id_type;
    typedef typename polygon_type::face_descripter   face_descripter;
    typedef typename polygon_type::edge_descripter   edge_descripter;
    typedef typename polygon_type::vertex_descripter vertex_descripter;
    typedef typename polygon_type::local_index_type  local_index_type;
    typedef typename polygon_type::barycentric_type  barycentric_type;

    typedef StructuralContainer<ParticleID, face_id_type, traits_type>
        structural_container_type;

protected:

    typedef typename ecell4::utils::pair_first_element_converter<
        Species, Species::serial_type, particle_id_set> serial_species_adapter;

public:

    ParticleSpacePolygonImpl(const boost::shared_ptr<polygon_type>& polygon)
        : base_type(), polygon_(polygon), container_(polygon->num_faces())
    {}
    virtual ~ParticleSpacePolygonImpl(){}

    virtual Integer num_species() const
    {
        return particle_pool_.size();
    }
    virtual bool has_particle(const Species& sp) const
    {
        return particle_pool_.count(sp.serial()) == 1;
    }

    virtual std::vector<Species> list_species() const
    {
        std::vector<Species> retval(particle_pool_.size());
        std::copy(boost::make_transform_iterator(
                particle_pool_.begin(), serial_species_adapter()),
            boost::make_transform_iterator(
                particle_pool_.end(), serial_species_adapter()),
            retval.begin());
        return retval;
    }

    const particle_container_type& particles() const
    {
        return particles_;
    }

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const;
    bool update_particle(const ParticleID& pid, const Particle& p);
    bool has_particle(const ParticleID& pid) const;
    void remove_particle(const ParticleID& pid);

    // polygon features

    face_id_type face_on(const ParticleID& pid) const;
    std::size_t  count_particle(const face_id_type& fid) const;
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type& fid);
    std::vector<std::pair<ParticleID, Particle> >
        get_particles(const face_id_type& fid) const;


    Integer num_particles() const {return particles_.size();}
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

    // for 3D
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const;

    // for 2D
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

protected:

    particle_container_type::iterator       find(const ParticleID&);
    particle_container_type::const_iterator find(const ParticleID&) const;

private:

    // and cell list?
    boost::shared_ptr<polygon_type> polygon_;
    particle_container_type         particles_;
    structural_container_type          container_;
    pid_to_idx_map_type             pid_to_idx_map_;
    per_species_particle_id_set     particle_pool_;
};

template<typename T_traits>
Integer
ParticleSpacePolygonImpl<T_traits>::num_particles(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    serial_species_adapter sadpt;

    for(per_species_particle_id_set::const_iterator
            iter(particle_pool_.begin()), end(particle_pool_.end());
            iter != end; ++iter)
    {
        if(sexp.match(sadpt(*iter))) retval += iter->second.size();
    }
    return retval;
}

template<typename T_traits>
Integer
ParticleSpacePolygonImpl<T_traits>::num_particles_exact(const Species& sp) const
{
    per_species_particle_id_set::const_iterator
        i(particle_pool_.find(sp.serial()));

    if (i == particle_pool_.end())
    {
        return 0;
    }
    return i->second.size();
}

template<typename T_traits>
Integer ParticleSpacePolygonImpl<T_traits>::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    serial_species_adapter   sadpt;
    for (per_species_particle_id_set::const_iterator i(particle_pool_.begin());
        i != particle_pool_.end(); ++i)
    {
        retval += sexp.count(sadpt(*i)) * i->second.size();
    }
    return retval;

}

template<typename T_traits>
Integer
ParticleSpacePolygonImpl<T_traits>::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

template<typename T_traits>
inline std::vector<std::pair<ParticleID, Particle> >
ParticleSpacePolygonImpl<T_traits>::list_particles() const
{
    return particles_;
}

template<typename T_traits>
std::vector<std::pair<ParticleID, Particle> >
ParticleSpacePolygonImpl<T_traits>::list_particles(const Species& sp) const
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

template<typename T_traits>
std::vector<std::pair<ParticleID, Particle> >
ParticleSpacePolygonImpl<T_traits>::list_particles_exact(const Species& sp) const
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

template<typename T_traits>
inline typename
ParticleSpacePolygonImpl<T_traits>::particle_container_type::iterator
ParticleSpacePolygonImpl<T_traits>::find(const ParticleID& k)
{
    const pid_to_idx_map_type::const_iterator p(pid_to_idx_map_.find(k));
    if(pid_to_idx_map_.end() == p)
        return particles_.end();
    return particles_.begin() + p->second;
}

template<typename T_traits>
inline typename
ParticleSpacePolygonImpl<T_traits>::particle_container_type::const_iterator
ParticleSpacePolygonImpl<T_traits>::find(const ParticleID& k) const
{
    const pid_to_idx_map_type::const_iterator p(pid_to_idx_map_.find(k));
    if(pid_to_idx_map_.end() == p)
        return particles_.end();
    return particles_.begin() + p->second;
}

template<typename T_traits>
std::pair<ParticleID, Particle>
ParticleSpacePolygonImpl<T_traits>::get_particle(const ParticleID& pid) const
{
    particle_container_type::const_iterator iter(find(pid));
    if(iter == particles_.end()) throw NotFound("particle not found");
    return *iter;
}

template<typename T_traits>
bool ParticleSpacePolygonImpl<T_traits>::update_particle(
        const ParticleID& pid, const Particle& p)
{
    if(container_.have(pid)) // 2D->3D
        container_.remove(pid);

    particle_container_type::iterator iter(find(pid));
    if(iter == particles_.end())
    {
        const std::size_t idx(particles_.size());
        pid_to_idx_map_[pid] = idx;
        particles_.push_back(std::make_pair(pid, p));
        return true;
    }
    else
    {
        assert(iter->first == pid);
        iter->second = p;
        return false;
    }
}

template<typename T_traits>
inline bool
ParticleSpacePolygonImpl<T_traits>::has_particle(const ParticleID& pid) const
{
    return find(pid) != particles_.end();
}

template<typename T_traits>
void ParticleSpacePolygonImpl<T_traits>::remove_particle(const ParticleID& pid)
{
    typename pid_to_idx_map_type::iterator iter(pid_to_idx_map_.find(pid));
    if(iter == pid_to_idx_map_.end())
        throw NotFound("particle not found");

    const std::size_t idx(iter->second), last_idx(particles_.size() - 1);
    if(idx != last_idx)
    {
        const std::pair<ParticleID, Particle>& last_elem(particles_.back());
        particles_[idx] = last_elem;
        pid_to_idx_map_[last_elem.first] = idx;
    }
    particles_.pop_back();
    pid_to_idx_map_.erase(pid);

    if(container_.have(pid))
    {
        container_.remove(pid);
    }
    return;
}

template<typename T_traits>
inline typename ParticleSpacePolygonImpl<T_traits>::face_id_type
ParticleSpacePolygonImpl<T_traits>::face_on(const ParticleID& pid) const
{
    if(container_.have(pid))
        return container_.structure_on(pid);
    throw NotFound("particle is not on a face");
}

template<typename T_traits>
inline std::size_t ParticleSpacePolygonImpl<T_traits>::count_particle(
        const face_id_type& fid) const
{
    return container_.elements_over(fid).size();
}

template<typename T_traits>
bool ParticleSpacePolygonImpl<T_traits>::update_particle(
        const ParticleID& pid, const Particle& p, const face_id_type& fid)
{
    particle_container_type::iterator iter(find(pid));
    if(iter == particles_.end()) // new particle
    {
        const std::size_t idx(particles_.size());
        pid_to_idx_map_[pid] = idx;
        particles_.push_back(std::make_pair(pid, p));
        container_.emplace(pid, fid);
        return true;
    }
    else if(container_.have(pid)) // 2D -> 2D
    {
        assert(iter->first == pid);
        iter->second = p;
        container_.update(pid, fid);
        return false;
    }
    else // 3D -> 2D
    {
        assert(iter->first == pid);
        iter->second = p;
        container_.emplace(pid, fid);
        return false;
    }
}

template<typename T_traits>
std::vector<std::pair<ParticleID, Particle> >
ParticleSpacePolygonImpl<T_traits>::get_particles(const face_id_type& fid) const
{
    const std::vector<ParticleID>& ids = container_.elements_over(fid);
    std::vector<std::pair<ParticleID, Particle> > retval(ids.size());
    for(std::size_t i=0; i<ids.size(); ++i)
    {
        retval[i] = this->get_particle(ids.at(i));
    }
    return retval;
}


template<typename T_traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpacePolygonImpl<T_traits>::list_particles_within_radius(
        const Real3& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > list;
    for(particle_container_type::const_iterator
        iter(particles_.begin()), end(particles_.end()); iter != end; ++iter)
    {
        const Real dist =
            length(pos - iter->second.position()) - iter->second.radius();
        if(dist < radius)
        {
            list.push_back(std::make_pair(*iter, dist));
        }
    }
    std::sort(list.begin(), list.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return list;
}

template<typename T_traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpacePolygonImpl<T_traits>::list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > list;
    for(particle_container_type::const_iterator
        iter(particles_.begin()), end(particles_.end()); iter != end; ++iter)
    {
        if(iter->first == ignore) continue;

        const Real dist =
            length(pos - iter->second.position()) - iter->second.radius();
        if(dist < radius)
        {
            list.push_back(std::make_pair(*iter, dist));
        }
    }
    std::sort(list.begin(), list.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return list;
}

template<typename T_traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpacePolygonImpl<T_traits>::list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > list;
    for(particle_container_type::const_iterator
        iter(particles_.begin()), end(particles_.end()); iter != end; ++iter)
    {
        if(iter->first == ignore1 || iter->first == ignore2) continue;

        const Real dist =
            length(pos - iter->second.position()) - iter->second.radius();
        if(dist < radius)
        {
            list.push_back(std::make_pair(*iter, dist));
        }
    }
    std::sort(list.begin(), list.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return list;
}

// for 2D
template<typename T_traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpacePolygonImpl<T_traits>::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    {// same face
        const std::vector<ParticleID>& ids = container_.elements_over(pos.second);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const particle_container_type::const_iterator found = find(*i);
            const Real dist = length(pos.first - found->second.position()) -
                              found->second.radius();
            if(dist < radius) retval.push_back(std::make_pair(*found, dist));
        }
    }

    boost::array<face_id_type, 3> const& adjacents =
        polygon_->adjacent_faces(pos.second);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::vector<ParticleID>& ids =
            container_.elements_over(adjacents[i]);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const particle_container_type::const_iterator found = find(*i);
            const Real dist = length(pos.first - found->second.position()) -
                              found->second.radius();
            if(dist < radius) retval.push_back(std::make_pair(*found, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

template<typename T_traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpacePolygonImpl<T_traits>::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    {// same face
        const std::vector<ParticleID>& ids = container_.elements_over(pos.second);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) continue;
            const particle_container_type::const_iterator found = find(*i);
            const Real dist = length(pos.first - found->second.position()) -
                              found->second.radius();
            if(dist < radius) retval.push_back(std::make_pair(*found, dist));
        }
    }

    boost::array<face_id_type, 3> const& adjacents =
        polygon_->adjacent_faces(pos.second);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::vector<ParticleID>& ids =
            container_.elements_over(adjacents[i]);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) continue;
            const particle_container_type::const_iterator found = find(*i);
            const Real dist = length(pos.first - found->second.position()) -
                              found->second.radius();
            if(dist < radius) retval.push_back(std::make_pair(*found, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

template<typename T_traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpacePolygonImpl<T_traits>::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    {// same face
        const std::vector<ParticleID>& ids = container_.elements_over(pos.second);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) continue;
            const particle_container_type::const_iterator found = find(*i);
            const Real dist = length(pos.first - found->second.position()) -
                              found->second.radius();
            if(dist < radius) retval.push_back(std::make_pair(*found, dist));
        }
    }

    boost::array<face_id_type, 3> const& adjacents =
        polygon_->adjacent_faces(pos.second);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::vector<ParticleID>& ids =
            container_.elements_over(adjacents[i]);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) continue;
            const particle_container_type::const_iterator found = find(*i);
            const Real dist = length(pos.first - found->second.position()) -
                              found->second.radius();
            if(dist < radius) retval.push_back(std::make_pair(*found, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

}// sgfrd
}// ecell4
#endif// ECELL4_SGFRD_PARTICLE_SPACE_POLYGON_IMPL
