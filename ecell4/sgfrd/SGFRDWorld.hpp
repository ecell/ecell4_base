#ifndef ECELL4_SGFRD_WORLD
#define ECELL4_SGFRD_WORLD
#include <ecell4/sgfrd/polygon_traits.hpp>
#include <ecell4/sgfrd/StructureRegistrator.hpp>
#include <ecell4/core/ParticleSpaceCellListImpl.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/Model.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/bind.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename T_polygon_traits>
class SGFRDWorld : public ecell4::Space
{
  public:
    typedef T_polygon_traits traits_type;
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

    typedef ecell4::Model model_type;

    typedef ParticleSpaceCellListImpl default_particle_space_type;
    typedef ParticleSpace particle_space_type;
    typedef typename particle_space_type::particle_container_type
        particle_container_type;
    typedef ecell4::SerialIDGenerator<ParticleID> particle_id_generator_type;
    typedef StructureRegistrator<ParticleID, face_id_type, traits_type>
        structure_registrator_type;

  public:

    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               const boost::shared_ptr<polygon_type>& polygon)
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          polygon_(polygon), registrator_(*polygon)
    {
        setup_descriptors(*polygon_);

        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        rng_->seed();
    }

    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               const boost::shared_ptr<polygon_type>& polygon,
               boost::shared_ptr<RandomNumberGenerator> rng)
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          rng_(rng), polygon_(polygon), registrator_(*polygon)
    {
        setup_descriptors(*polygon_);
    }

    ~SGFRDWorld(){}

    boost::shared_ptr<RandomNumberGenerator> const& rng() {return this->rng_;}
    boost::shared_ptr<polygon_type> const& polygon() const {return polygon_;}

    const Real t()                  const {return ps_->t();}
    void       set_t(const Real& t)       {return ps_->set_t(t);}
    const Real volume()             const {return ps_->volume();}
    Real get_value(const Species& sp)       const {return ps_->get_value(sp);}
    Real get_value_exact(const Species& sp) const {return ps_->get_value_exact(sp);}

    Integer num_species() const {return ps_->num_species();}
    bool has_species(const Species& sp) const {return ps_->has_species(sp);}
    std::vector<Species> list_species() const {return ps_->list_species();}

// CellListImpl stuff
//     Real3 const& edge_lengths() const {return ps_->edge_lengths();}
//     Real3 const& cell_sizes()   const {return ps_->cell_sizes();}
//     Integer3     matrix_sizes() const {return ps_->matrix_sizes();}
//     void reset(const Real3& edge_lengths) {ps_->reset(edge_lengths);}

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p);
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p, const face_id_type& fid);

    bool update_particle(const ParticleID& pid, const Particle& p);
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type& fid);

    // this also removes particle if it is on surface
    void
    remove_particle(const ParticleID& pid);
    void
    remove_particle(const ParticleID& pid, const face_id_type& fid);

    std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const {return ps_->get_particle(pid);}
    bool
    has_particle(const ParticleID& pid) const {return ps_->has_particle(pid);}
    bool
    has_particle(const face_id_type& fid) const;
    Integer
    num_particle(const face_id_type& fid) const;

    Integer
    num_particles()                        const {return ps_->num_particles();}
    Integer
    num_particles(const Species& sp)       const {return ps_->num_particles(sp);}
    Integer
    num_particles_exact(const Species& sp) const {return ps_->num_particles_exact(sp);}
    Integer
    num_molecules(const Species& sp)       const {return ps_->num_molecules(sp);}
    Integer
    num_molecules_exact(const Species& sp) const {return ps_->num_molecules_exact(sp);}

    std::vector<std::pair<ParticleID, Particle> >
    list_particles() const {return ps_->list_particles();}
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const {return ps_->list_particles(sp);}
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const {return ps_->list_particles_exact(sp);}

    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const face_id_type& fid) const;
    std::vector<ParticleID> const&
    list_particleIDs(const face_id_type& fid) const;

    face_id_type get_face_id(const ParticleID& pid) const;

    void save(const std::string& fname) const
    {
        throw NotImplemented("SGFRDWorld::save");
    }
    void load(const std::string& fname) const
    {
        throw NotImplemented("SGFRDWorld::load");
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        throw NotImplemented("SGFRDWorld::save_hdf5");
    }

    void load_hdf5(const H5::Group& root)
    {
        throw NotImplemented("SGFRDWorld::load_hdf5");
    }
#endif

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

    particle_container_type const& particles() const {return ps_->particles();}


    void bind_to(boost::shared_ptr<model_type> model)
    {
        if (boost::shared_ptr<model_type> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to BDWorld"
                    << std::endl;
            }
        }
        model_ = model;
    }

    boost::shared_ptr<model_type> lock_model() const
    {
        return model_.lock();
    }

  private:

    boost::scoped_ptr<particle_space_type>   ps_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
    boost::weak_ptr<Model>                   model_;
    boost::shared_ptr<polygon_type>          polygon_;
    structure_registrator_type               registrator_;
    particle_id_generator_type               pidgen_;
};

template<typename traits>
std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld<traits>::new_particle(const Particle& p)
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

template<typename traits>
std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld<traits>::new_particle(const Particle& p, const face_id_type& fid)
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

template<typename traits>
inline bool
SGFRDWorld<traits>::update_particle(const ParticleID& pid, const Particle& p)
{
    return ps_->update_particle(pid, p);
}

template<typename traits>
inline bool
SGFRDWorld<traits>::update_particle(const ParticleID& pid, const Particle& p,
                                    const face_id_type& fid)
{
    if(registrator_.have(pid))
    {
        registrator_.update(pid, fid);
    }
    else
    {
        registrator_.emplace(pid, fid);
    }
    return ps_->update_particle(pid, p);
}

template<typename traits>
inline void
SGFRDWorld<traits>::remove_particle(const ParticleID& pid)
{
    if(registrator_.have(pid)) registrator_.remove(pid);
    return ps_->remove_particle(pid);
}

template<typename traits>
inline void
SGFRDWorld<traits>::remove_particle(
        const ParticleID& pid, const face_id_type& fid)
{
    registrator_.remove(pid, fid);
    return ps_->remove_particle(pid);
}

template<typename traits>
inline bool
SGFRDWorld<traits>::has_particle(const face_id_type& fid) const
{
    return registrator_.elements_over(fid).size() > 0;
}

template<typename traits>
inline Integer
SGFRDWorld<traits>::num_particle(const face_id_type& fid) const
{
    return registrator_.elements_over(fid).size();
}

template<typename traits>
std::vector<std::pair<ParticleID, Particle> >
SGFRDWorld<traits>::list_particles(const face_id_type& fid) const
{
    std::vector<ParticleID> const& pids = registrator_.elements_over(fid);
    std::vector<std::pair<ParticleID, Particle> > retval(pids.size());
    std::transform(pids.begin(), pids.end(), retval.begin(),
            boost::bind(&SGFRDWorld<traits>::get_particle, this, _1));
    return retval;
}

template<typename traits>
inline std::vector<ParticleID> const&
SGFRDWorld<traits>::list_particleIDs(const face_id_type& fid) const
{
    return registrator_.elements_over(fid);
}

template<typename traits>
inline typename SGFRDWorld<traits>::face_id_type
SGFRDWorld<traits>::get_face_id(const ParticleID& pid) const
{
    return registrator_.structure_on(pid);
}


template<typename traits>
inline std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld<traits>::list_particles_within_radius(
            const Real3& pos, const Real& radius) const
{
    return ps_->list_particles_within_radius(pos, radius);
}

template<typename traits>
inline std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld<traits>::list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore) const
{
    return ps_->list_particles_within_radius(pos, radius, ignore);
}

template<typename traits>
inline std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld<traits>::list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
{
    return ps_->list_particles_within_radius(pos, radius, ignore1, ignore2);
}

template<typename traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld<traits>::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(typename std::vector<ParticleID>::const_iterator
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
    for(typename std::vector<face_id_type>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = registrator_.elements_over(*iter);
        for(typename std::vector<ParticleID>::const_iterator
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

template<typename traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld<traits>::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(typename std::vector<ParticleID>::const_iterator
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
    for(typename std::vector<face_id_type>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(typename std::vector<ParticleID>::const_iterator
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

template<typename traits>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld<traits>::list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(typename std::vector<ParticleID>::const_iterator
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
    for(typename std::vector<face_id_type>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(typename std::vector<ParticleID>::const_iterator
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
#endif // ECELL4_SGFRD_WORLD
