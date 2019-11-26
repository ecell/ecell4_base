#ifndef ECELL4_SGFRD_WORLD
#define ECELL4_SGFRD_WORLD
#include "StructureRegistrator.hpp"
#include "ReactionInfo.hpp"
#include <ecell4/core/extras.hpp>
#include <ecell4/core/WorldInterface.hpp>
#include <ecell4/core/ParticleSpaceCellListImpl.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/STLFileIO.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/Segment.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/collision.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared.hpp>
#include <boost/weak_ptr.hpp>
#include <array>

namespace ecell4
{
namespace sgfrd
{

struct MoleculeInfo
{
    const Real radius;
    const Real D;
};

class SGFRDWorld
    : public ecell4::WorldInterface
{
  public:
    typedef ecell4::Polygon  polygon_type;
    typedef polygon_type::FaceID   FaceID;
    typedef polygon_type::EdgeID   EdgeID;
    typedef polygon_type::VertexID VertexID;
    typedef Barycentric barycentric_type;

    typedef ecell4::sgfrd::MoleculeInfo molecule_info_type;
    typedef ecell4::Model model_type;

    typedef ParticleSpaceCellListImpl default_particle_space_type;
    typedef ParticleSpace particle_space_type;
    typedef particle_space_type::particle_container_type
        particle_container_type;
    typedef ecell4::SerialIDGenerator<ParticleID> particle_id_generator_type;
    typedef StructureRegistrator<ParticleID, FaceID>
        structure_registrator_type;

  public:

    // !rng && !polygon_file
    SGFRDWorld(const Real3&    edge_lengths =    Real3(1, 1, 1),
               const Integer3& matrix_sizes = Integer3(3, 3, 3))
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          polygon_(boost::make_shared<Polygon>(edge_lengths, matrix_sizes)),
          registrator_(*polygon_)
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        rng_->seed();

        this->prepare_restrictions();
    }

    // rng && !polygon_file
    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               boost::shared_ptr<RandomNumberGenerator> rng)
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          rng_(rng),
          polygon_(boost::make_shared<Polygon>(edge_lengths, matrix_sizes)),
          registrator_(*polygon_)
    {
        this->prepare_restrictions();
    }

    // !rng && polygon_file
    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               const std::string& polygon_file, const STLFormat fmt)
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          polygon_(boost::make_shared<Polygon>(
                      read_polygon(polygon_file, fmt, edge_lengths))),
          registrator_(*polygon_)
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        rng_->seed();

        this->prepare_restrictions();
    }

    // rng && polygon_file
    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               boost::shared_ptr<RandomNumberGenerator> rng,
               const std::string& polygon_file, const STLFormat fmt)
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          rng_(rng), polygon_(boost::make_shared<Polygon>(
                      read_polygon(polygon_file, fmt, edge_lengths))),
          registrator_(*polygon_)
    {
        this->prepare_restrictions();
    }

    // !rng && polygon
    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               const boost::shared_ptr<Polygon>& poly)
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          polygon_(poly), registrator_(*polygon_)
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        rng_->seed();

        this->prepare_restrictions();

        write_polygon("tmp.stl", STLFormat::Ascii, *polygon());
    }

    // rng && polygon
    SGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               boost::shared_ptr<RandomNumberGenerator> rng,
               const boost::shared_ptr<Polygon>& poly)
        : ps_(new default_particle_space_type(edge_lengths, matrix_sizes)),
          rng_(rng), polygon_(poly), registrator_(*polygon_)
    {
        this->prepare_restrictions();
    }

    SGFRDWorld(const std::string& filename) // from HDF5
        : ps_(new default_particle_space_type(Real3(1, 1, 1))),
          polygon_(boost::make_shared<Polygon>(Real3(1, 1, 1))),
          registrator_(*polygon_)
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }


    ~SGFRDWorld() override = default;

    boost::shared_ptr<RandomNumberGenerator> const& rng() const noexcept {return this->rng_;}
    boost::shared_ptr<RandomNumberGenerator>&       rng()       noexcept {return this->rng_;}
    boost::shared_ptr<polygon_type> const& polygon() const {return polygon_;}

    // -----------------------------------------------------------------------
    // WorldInterface

    const Real t() const override
    {
        return ps_->t();
    }
    void set_t(const Real& t) override
    {
        return ps_->set_t(t);
    }

    void save(const std::string& filename) const override
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        pidgen_.save(fout.get());

        boost::scoped_ptr<H5::Group>
            group1(new H5::Group(fout->createGroup("ParticleSpace")));
        ps_->save_hdf5(group1.get());

        boost::scoped_ptr<H5::Group>
            group2(new H5::Group(fout->createGroup("Polygon")));
        this->polygon_->save_hdf5(group2.get());

        extras::save_version_information(
            fout.get(), std::string("ecell4-sgfrd-") + std::string(VERSION_INFO));
        return;
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
    }

    void load(const std::string& filename) override
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));

        const std::string required = "ecell4-sgfrd-0.0";
        try
        {
            const std::string version = extras::load_version_information(*fin);
            if (!extras::check_version_information(version, required))
            {
                std::stringstream ss;
                ss << "The version of the given file [" << version
                    << "] is too old. [" << required << "] or later is required.";
                throw NotSupported(ss.str());
            }
        }
        catch(H5::GroupIException not_found_error)
        {
            throw NotFound("No version information was found.");
        }

        const H5::Group group1(fin->openGroup("ParticleSpace"));
        ps_->load_hdf5(group1);

        const H5::Group group2(fin->openGroup("Polygon"));
        polygon_->load_hdf5(group2);

        pidgen_.load(*fin);
        rng_->load(*fin);

        // restore polygon-particle relationships.
        // --------------------------------------------------------------------
        // The above code reads 3D positions of particles. But the information
        // about which particle is on which face are not restored.
        // The nice thing is that all the restored particles locate precisely
        // (within the limit of numerical error) on the faces unless the polygon
        // shape was changed by restoring.

        for(auto pidp : this->list_particles())
        {
            const auto& fp = this->find_face(pidp.second.position());
            if(!fp)
            {
                throw std::invalid_argument("Particle does not locate on any face");
            }
            pidp.second.position() = fp->first;
            this->update_particle(pidp.first, pidp.second, fp->second);
        }
        return;
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
    }

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    MoleculeInfo get_molecule_info(const Species& sp) const
    {
        Real radius(0.0), D(0.0);

        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            auto const sp_ = bound_model->apply_species_attributes(sp);
            if (sp_.has_attribute("radius"))
            {
                radius = sp_.get_attribute_as<Real>("radius");
            }
            if (sp_.has_attribute("D"))
            {
                D = sp_.get_attribute_as<Real>("D");
            }
        }

        {
            if (sp.has_attribute("radius"))
            {
                radius = sp.get_attribute_as<Real>("radius");
            }
            if (sp.has_attribute("D"))
            {
                D = sp.get_attribute_as<Real>("D");
            }
        }

        if (radius <= 0.0)
        {
            std::stringstream msg;
            msg << "A particle with invalid size [" << radius << "] was given.";
            throw IllegalArgument(msg.str());
        }

        MoleculeInfo info = {radius, D};
        return info;
    }

    const Real volume() const override
    {
        return ps_->volume();
    }

    bool has_species(const Species& sp) const override
    {
        for(auto&& item : ps_->list_species())
        {
            if(item == sp){return true;}
        }
        return false;
    }

    std::vector<Species> list_species() const
    {
        return ps_->list_species();
    }

    Integer num_molecules(const Species& sp) const override
    {
        return ps_->num_molecules(sp);
    }

    Integer num_molecules_exact(const Species& sp) const override
    {
        return ps_->num_molecules_exact(sp);
    }

    Real get_value(const Species& sp)       const override
    {
        return ps_->get_value(sp);
    }
    Real get_value_exact(const Species& sp) const override
    {
        return ps_->get_value_exact(sp);
    }

    const Real3& edge_lengths() const override
    {
        return ps_->edge_lengths();
    }

    Integer num_particles() const override
    {
        return ps_->num_particles();
    }
    Integer num_particles(const Species& sp) const override
    {
        return ps_->num_particles(sp);
    }
    Integer num_particles_exact(const Species& sp) const override
    {
        return ps_->num_particles_exact(sp);
    }

    bool has_particle(const ParticleID& pid) const override
    {
        return ps_->has_particle(pid);
    }
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const override
    {
        return ps_->get_particle(pid);
    }

    std::pair<FaceID, Barycentric>
    get_surface_position(const ParticleID& pid) const
    {
        const auto  fid = this->registrator_.structure_on(pid);
        const auto& tri = this->polygon_->triangle_at(fid);
        const auto  pos = this->ps_->get_particle(pid).second.position();
        return std::make_pair(fid, to_barycentric(pos, tri));
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles() const override
    {
        return ps_->list_particles();
    }
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const override
    {
        return ps_->list_particles(sp);
    }
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const override
    {
        return ps_->list_particles_exact(sp);
    }

    std::vector<std::pair<ParticleID, std::pair<FaceID, Barycentric>>>
    list_surface_positions() const
    {
        std::vector<std::pair<ParticleID, std::pair<FaceID, Barycentric>>> v;
        v.reserve(this->num_particles());
        for(auto&& pp : this->list_particles())
        {
            v.emplace_back(pp.first, this->get_surface_position(pp.first));
        }
        return v;
    }
    std::vector<std::pair<ParticleID, std::pair<FaceID, Barycentric>>>
    list_surface_positions(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, std::pair<FaceID, Barycentric>>> v;
        v.reserve(this->num_particles(sp));
        for(auto&& pp : this->list_particles(sp))
        {
            v.emplace_back(pp.first, this->get_surface_position(pp.first));
        }
        return v;
    }
    std::vector<std::pair<ParticleID, std::pair<FaceID, Barycentric>>>
    list_surface_positions_exact(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, std::pair<FaceID, Barycentric>>> v;
        v.reserve(this->num_particles_exact(sp));
        for(auto&& pp : this->list_particles_exact(sp))
        {
            v.emplace_back(pp.first, this->get_surface_position(pp.first));
        }
        return v;
    }

    // -----------------------------------------------------------------------
    // ParticleSpaceInterface

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p);
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p, const FaceID& fid);

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Species& sp, const Real3& pos)
    {
        const auto info = this->get_molecule_info(sp);
        return this->new_particle(Particle(sp, pos, info.radius, info.D));
    }
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Species& sp, const FaceID& fid, const Barycentric& bary)
    {
        const auto info = this->get_molecule_info(sp);

        const auto& tri = this->polygon_->triangle_at(fid);
        const auto  pos = to_absolute(bary, tri);

        return this->new_particle(Particle(sp, pos, info.radius, info.D), fid);
    }
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Species& sp, const std::pair<FaceID, Barycentric>& sfp)
    {
        return this->new_particle(sp, sfp.first, sfp.second);
    }

    std::pair<std::pair<ParticleID, Particle>, bool>
    throw_in_particle(const Species& sp);

    void add_molecules(const Species& sp, const Integer& num);
    void add_molecules(const Species& sp, const Integer& num,
                       const boost::shared_ptr<Shape> shape);

    void remove_molecules(const Species& sp, const Integer& num)
    {
        if (num < 0)
        {
            throw std::invalid_argument(
                "The number of molecules must be positive.");
        }

        auto particles(list_particles(sp));
        const Integer num_particles(particles.size());
        if (num_particles < num)
        {
            throw std::invalid_argument(
                "The number of molecules cannot be negative.");
        }

        shuffle((*rng_), particles);
        for (auto i(particles.begin()); i != particles.begin() + num; ++i)
        {
            remove_particle((*i).first);
        }
        return;
    }

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        // Note: condition in `if` statements can have a declarator.
        //       cf. N3337 section 6.4 "Selection statements"
        if(const auto pfid = this->find_face(p.position()))
        {
            Particle p_(p);
            p_.position() = pfid->first;

            return this->update_particle(pid, p_, pfid->second);
        }
        throw std::invalid_argument("[error] SGFRDWorld::update_particle: "
                "particle locates distant from polygon");
    }
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const FaceID& fid)
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

    // this also removes particle if it is on surface
    void remove_particle(const ParticleID& pid)
    {
        if(registrator_.have(pid))
        {
            registrator_.remove(pid);
        }
        return ps_->remove_particle(pid);
    }
    void remove_particle(const ParticleID& pid, const FaceID& fid)
    {
        registrator_.remove(pid, fid);
        return ps_->remove_particle(pid);
    }
    bool
    has_particle(const FaceID& fid) const
    {
        return registrator_.elements_over(fid).size() > 0;
    }
    Integer
    num_particle(const FaceID& fid) const
    {
        return registrator_.elements_over(fid).size();
    }

    Triangle get_triangle(const ParticleID& pid) const
    {
        return this->polygon_->triangle_at(this->registrator_.structure_on(pid));
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const FaceID& fid) const
    {
        const std::vector<ParticleID>& pids = registrator_.elements_over(fid);
        std::vector<std::pair<ParticleID, Particle>> retval(pids.size());

        std::transform(pids.begin(), pids.end(), retval.begin(),
            [this](const ParticleID& pid) -> std::pair<ParticleID, Particle> {
                return this->get_particle(pid);
            });
        return retval;
    }
    std::vector<ParticleID> const&
    list_particleIDs(const FaceID& fid) const
    {
        return registrator_.elements_over(fid);
    }

    bool is_on_face(const ParticleID& pid) const
    {
        return registrator_.have(pid);
    }

    FaceID get_face_id(const ParticleID& pid) const
    {
        return registrator_.structure_on(pid);
    }

    Real distance_sq(const Real3& lhs, const Real3& rhs)
    {
        const auto pf1 = this->find_face(lhs);
        if(pf1)
        {
            const auto pf2 = this->find_face(rhs);
            if(pf2)
            {
                return this->distance_sq(*pf1, *pf2);
            }
        }
        throw std::invalid_argument("[error] SGFRDWorld::distance_sq: "
                "particle locates distant from polygon");
    }
    Real distance(const Real3& lhs, const Real3& rhs)
    {
        const auto pf1 = this->find_face(lhs);
        if(pf1)
        {
            const auto pf2 = this->find_face(rhs);
            if(pf2)
            {
                return this->distance(*pf1, *pf2);
            }
        }
        throw std::invalid_argument("[error] SGFRDWorld::distance: "
                "particle locates distant from polygon");
    }

    template<typename str1T, typename str2T>
    Real distance_sq(const std::pair<Real3, str1T>& lhs,
                     const std::pair<Real3, str2T>& rhs)
    {
        return ecell4::polygon::distance_sq(*polygon_, lhs, rhs);
    }
    template<typename str1T, typename str2T>
    Real distance(const std::pair<Real3, str1T>& lhs,
                  const std::pair<Real3, str2T>& rhs)
    {
        return ecell4::polygon::distance(*polygon_, lhs, rhs);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius) const
    {
        if(const auto pf = this->find_face(pos))
        {
            return this->list_particles_within_radius(*pf, radius);
        }
        throw std::invalid_argument("[error] "
                "SGFRDWorld::list_particles_within_radius: "
                "particle locates distant from polygon");
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore) const
    {
        if(const auto pf = this->find_face(pos))
        {
            return this->list_particles_within_radius(*pf, radius, ignore);
        }
        throw std::invalid_argument("[error] "
                "SGFRDWorld::list_particles_within_radius: "
                "particle locates distant from polygon");
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        if(const auto pf = this->find_face(pos))
        {
            return this->list_particles_within_radius(
                    *pf, radius, ignore1, ignore2);
        }
        throw std::invalid_argument("[error] "
                "SGFRDWorld::list_particles_within_radius: "
                "particle locates distant from polygon");
    }

    // for 2D
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const std::pair<Real3, FaceID>& pos, const Real& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const std::pair<Real3, FaceID>& pos, const Real& radius,
            const ParticleID& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const std::pair<Real3, FaceID>& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const;

    std::vector<std::pair<VertexID, Real> > list_vertices_within_radius(
            const std::pair<Real3, FaceID>& pos, const Real& radius) const
    {
        const Real threshold_sq = radius * radius;
        std::vector<std::pair<VertexID, Real> > retval;
        const std::vector<VertexID> candidates =
            this->polygon_->neighbor_vertices_of(pos.second);
        for(std::vector<VertexID>::const_iterator
                i(candidates.begin()), e(candidates.end()); i!=e; ++i)
        {
            const VertexID vid = *i;
            const Real dist_sq = ecell4::polygon::distance_sq(*(this->polygon_),
                    pos, std::make_pair(this->polygon_->position_at(vid), vid));
            if(dist_sq < threshold_sq)
            {
                retval.push_back(std::make_pair(vid, std::sqrt(dist_sq)));
            }
        }
        std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<VertexID, Real>());
        return retval;
    }

    // return false if overlap exists. for 3D. FIXME: speedup
    bool check_no_overlap(const Real3& pos, const Real& radius) const
    {
        return this->list_particles_within_radius(pos, radius).empty();
    }
    bool check_no_overlap(const Real3& pos, const Real& radius,
            const ParticleID& ignore) const
    {
        return this->list_particles_within_radius(pos, radius, ignore).empty();
    }
    bool check_no_overlap(const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return this->list_particles_within_radius(pos, radius, ignore1, ignore2
                ).empty();
    }

    // return false if overlap exists.
    bool check_no_overlap(
            const std::pair<Real3, FaceID>& pos, const Real& radius) const;
    bool check_no_overlap(
            const std::pair<Real3, FaceID>& pos, const Real& radius,
            const ParticleID& ignore) const;
    bool check_no_overlap(
            const std::pair<Real3, FaceID>& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const;

    particle_container_type const& particles() const {return ps_->particles();}

    Real3 periodic_transpose(const Real3& pos1, const Real3& pos2) const
    {
        return this->polygon_->periodic_transpose(pos1, pos2);
    }
    Real3 apply_boundary(const Real3& pos) const
    {
        return this->ps_->apply_boundary(pos);
    }

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

    std::array<ecell4::Segment, 6> const& barrier_at(const FaceID& fid) const
    {
        return this->barriers_.at(fid);
    }

    Real estimated_possible_largest_particle_radius() const noexcept
    {
        return this->estimated_possible_largest_particle_radius_;
    }

  private:

    // the tolerance is relative to edge_lengths.
    boost::optional<std::pair<Real3, FaceID>>
    find_face(const Real3& pos, const Real tolerance = 1e-3) const
    {
        const auto& width = this->edge_lengths();
        const auto tol = tolerance * std::min(width[0], std::min(width[1], width[2]));
        const auto tol2 = tol * tol;

        Real min_distance = std::numeric_limits<Real>::infinity();
        boost::optional<std::pair<Real3, FaceID>> nearest = boost::none;
        for(const auto& fid : this->polygon_->list_face_ids())
        {
            const auto& tri = this->polygon_->triangle_at(fid);
            const auto dist =
                ecell4::collision::distance_sq_point_triangle(pos, tri);
            if(dist <= tol2)
            {
                if(!nearest || dist < min_distance)
                {
                    min_distance = dist;
                    nearest      = std::make_pair(pos, fid);
                }
            }
        }
        return nearest;
    }

    inline static Real3 normalize(const Real3& v) throw()
    {
        return v * (1.0 / std::sqrt(length_sq(v)));
    }

    void prepare_restrictions()
    {
        // To avoid edge cases, it calculates the maximum size of particle.
        // Also, to avoid overlap between shells, it calculates a bisector of
        // each angle in triangle.
        Real min_altitude = std::numeric_limits<Real>::max();
        for(const auto& fid : polygon_->list_face_ids())
        {
            const auto& tri = polygon_->triangle_at(fid);

            // Estimate largest particle radius possible.
            const auto S = tri.area();
            min_altitude = std::min(min_altitude, 2.0 * S / tri.length_of_edge_at(0));
            min_altitude = std::min(min_altitude, 2.0 * S / tri.length_of_edge_at(1));
            min_altitude = std::min(min_altitude, 2.0 * S / tri.length_of_edge_at(2));

            // calculate boundary for shell size
            const auto& edges = polygon_->edges_of(fid);
            const auto& vtxs  = polygon_->vertices_of(fid);
            std::array<Segment, 6> segments;
            for(std::size_t i=0; i<3; ++i)
            {
                // vi1   ei1  vi0  |
                //     <-----.     |
                //     \    ^ \    |
                // ei2  \  /ei0\   |
                //       v/_____\  |
                //      vi2        |

                const auto  ei0 = polygon_->opposite_of(edges.at(i));
                const auto  ei1 = polygon_->next_of(ei0);
                const auto  ei2 = polygon_->next_of(ei1);
                const auto lei0 = polygon_->length_of(ei0);
                const auto lei1 = polygon_->length_of(ei1);
                const auto lei2 = polygon_->length_of(ei2);
                const auto dei1 = polygon_->direction_of(ei1);
                const auto dei2 = polygon_->direction_of(ei2);

                const auto  vi0 = polygon_->target_of(ei0);
                const auto  vi1 = polygon_->target_of(ei1);
                const auto  vi2 = polygon_->target_of(ei2);

                assert(vi0 == vtxs[i]);
                assert(vi2 == vtxs[i==2?0:i+1]);

                const auto pvi0 = tri.vertices()[i];
                const auto pvi1 = this->periodic_transpose(polygon_->position_at(vi1), pvi0);
                const auto pvi2 = tri.vertices()[(i==2)?0:i+1];

                const auto dst0 = pvi1 + dei2 * (lei1 / (lei1 + lei0));
                const auto dst2 = pvi0 + dei1 * (lei0 / (lei0 + lei2));

                segments[2*i  ] = Segment(dst0, pvi0);
                segments[2*i+1] = Segment(dst2, pvi2);
            }
            this->barriers_[fid] = segments;
        }
        this->estimated_possible_largest_particle_radius_ = min_altitude * 0.5;
        return;
    }

  private:

    boost::scoped_ptr<particle_space_type>   ps_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
    boost::weak_ptr<Model>                   model_;
    boost::shared_ptr<polygon_type>          polygon_;
    structure_registrator_type               registrator_;
    particle_id_generator_type               pidgen_;
    Real estimated_possible_largest_particle_radius_;

    // this contains the edges that correspond to the developed neighbor faces.
    //
    //        /\
    //     > /__\ <
    //      /\* /\
    //   > /__\/__\ < these edges
    //      ^    ^
    boost::container::flat_map<FaceID, std::array<ecell4::Segment, 6> > barriers_;
};


}// sgfrd
}// ecell4
#endif // ECELL4_SGFRD_WORLD
