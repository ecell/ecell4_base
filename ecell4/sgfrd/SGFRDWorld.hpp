#ifndef ECELL4_SGFRD_WORLD
#define ECELL4_SGFRD_WORLD
#include <ecell4/sgfrd/StructureRegistrator.hpp>
#include <ecell4/core/ParticleSpaceCellListImpl.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/Context.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename T_polygon_traits>
class SGFRDWorld : public ecell4::Space
{
  public:
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

    typedef ParticleSpaceCellListImpl particle_space_type;
    typedef StructureRegistrator<ParticleID, face_id_type, traits_type>
        structure_registrator_type;

  public:

    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               const polygon_type& polygon)
        : pcon_(new particle_space_type(edge_lengths, matrix_sizes)),
          polygon_(polygon), registrator_(polygon.num_faces())
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        rng_->seed();
    }

    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               const polygon_type& polygon,
               boost::shared_ptr<RandomNumberGenerator> rng)
        : pcon_(new particle_space_type(edge_lengths, matrix_sizes)),
          rng_(rng), polygon_(polygon), registrator_(polygon.num_faces())
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        rng_->seed();
    }

    ~SGFRDWorld(){}

    const Real t()                  const {return ps_->t();}
    void       set_t(const Real& t)       {return ps_->set_t(t);}
    const Real volume()             const {return ps_->volume();}
    Real get_value(const Species& sp)       const {return ps_->get_value(sp);}
    Real get_value_exact(const Species& sp) const {return ps_->get_value_exact(sp);}

    Integer num_species() const {return ps_->num_species();}
    bool has_species(const Species& sp) const {return ps_->has_species(sp);}
    std::vector<Species> list_species() const {return ps_->list_species();}

    Real3 const& edge_lengths() const {return ps_->edge_lengths();}
    Real3 const& cell_sizes()   const {return ps_->cell_sizes();}
    Integer3     matrix_sizes() const {return ps_->matrix_sizes();}

    void reset(const Real3& edge_lengths) {ps_->reset(edge_lengths);}
    bool update_particle(const ParticleID& pid, const Particle& p);
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type& fid);

    std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const {return ps_->get_particle(pid);}
    bool
    has_particle(const ParticleID& pid) const {return ps_->has_particle(pid);}
    void
    remove_particle(const ParticleID& pid) const;

    Integer num_particles() const {return particles_.size();}
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

    //TODO 2D accessors...

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

    particle_container const& particles() const {return ps_->particles();}

  private:

    boost::scoped_ptr<particle_space_type>   ps_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
    polygon_type                             polygon_;
    structure_registrator_type               registrator_;
};


}// sgfrd
}// ecell4
#endif // ECELL4_SGFRD_WORLD
