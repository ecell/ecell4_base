#ifndef ECELL4_NGFRD_NGFRD_WORLD_HPP
#define ECELL4_NGFRD_NGFRD_WORLD_HPP

#include <memory>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/ParticleSpaceRTreeImpl.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/WorldInterface.hpp>
#include <ecell4/ngfrd/PolygonContainer.hpp>

namespace ecell4
{
namespace ngfrd
{

struct MoleculeInfo
{
    const Real radius;
    const Real D;
};

class NGFRDWorld final
    : public WorldInterface
{
public:
    using molecule_info_type      = MoleculeInfo;
    using particle_space_type     = ParticleSpaceRTreeImpl;
    using particle_container_type = particle_space_type::particle_container_type;
    using polygon_container_type  = PolygonContainer<ParticleID>;

public:

    NGFRDWorld(const Real3&    edge_lengths = Real3(1, 1, 1),
               const Integer3& matrix_sizes = Integer3(3, 3, 3))
        : ps_(new particle_space_type(edge_lengths, matrix_sizes)),
          polygon_(std::make_shared<Polygon>(edge_lengths, matrix_sizes))
    {
        rng_ = std::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();

        this->prepare_restrictions();
    }

    NGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               std::shared_ptr<RandomNumberGenerator> rng,
               std::shared_ptr<Polygon>               poly)
        : ps_(new particle_space_type(edge_lengths, matrix_sizes)),
          rng_(std::move(rng)), polygon_(std::move(poly))
    {
        this->prepare_restrictions();
    }

    NGFRDWorld(const std::string& filename)
        : ps_(new particle_space_type(Real3(1, 1, 1)))
    {
        rng_ = std::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }

    ~NGFRDWorld() override = default;

    // ------------------------------------------------------------------------

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p);

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p, const FaceID& fid);

    bool update_particle(const ParticleID& pid, const Particle& p);
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const FaceID& fid);

    void remove_particle(const ParticleID&);
    bool has_particle(const ParticleID&);

    // ------------------------------------------------------------------------

    boost::optional<FaceID> on_which_face(const ParticleID& pid) const noexcept
    {
        return poly_con_.on_which_face(pid);
    }

    boost::optional<std::vector<ObjectID> const&>
    particles_on(const FaceID& fid) const noexcept
    {
        return poly_con_.objects_on(fid);
    }

    // ------------------------------------------------------------------------

    bool has_overlapping_particle(const Real3&) const;
    bool has_overlapping_particle(const Real3&, const ParticleID& ignore) const;
    bool has_overlapping_particle(const Real3&, const ParticleID& ignore1,
            const ParticleID& ignore2) const;

    bool has_overlapping_particle(const std::pair<Real3, FaceID>&) const;
    bool has_overlapping_particle(const std::pair<Real3, FaceID>&,
            const ParticleID&) const;
    bool has_overlapping_particle(const std::pair<Real3, FaceID>&,
            const ParticleID&, const ParticleID&) const;

    // ------------------------------------------------------------------------
    // list_particles_within_radius

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& center, const Real radius) const
    {
        return list_particles_within_radius_impl(center, radius,
                [](const ParticleID&) noexcept -> bool {return false;});
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& center, const Real radius,
            const ParticleID& ignore) const
    {
        return list_particles_within_radius_impl(center, radius,
                [&](const ParticleID& pid) noexcept -> bool {
                    return pid == ignore;
                });
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& center, const Real radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return list_particles_within_radius_impl(center, radius,
                [&](const ParticleID& pid) noexcept -> bool {
                    return pid == ignore1 || pid == ignore2;
                });
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(
            const std::pair<Real3, FaceID>& center, const Real radius) const
    {
        return list_particles_within_radius_impl(center, radius,
                [](const ParticleID&) noexcept -> bool {return false;});
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(
            const std::pair<Real3, FaceID>& center, const Real radius,
            const ParticleID& ignore) const
    {
        return list_particles_within_radius_impl(center, radius,
                [&](const ParticleID& pid) noexcept -> bool {
                    return pid == ignore;
                });
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(
            const std::pair<Real3, FaceID>& center, const Real radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return list_particles_within_radius_impl(center, radius,
                [&](const ParticleID& pid) noexcept -> bool {
                    return pid == ignore1 || pid == ignore2;
                });
    }

    Polygon&       polygon()       noexcept {return *polygon_;}
    Polygon const& polygon() const noexcept {return *polygon_;}

    PeriodicBoundary boundary() const noexcept
    {
        return PeriodicBoundary(edge_lengths);
    }

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
        throw NotImplemented("NGFRD::save(): TODO");
    }
    void load(const std::string& filename) override
    {
        throw NotImplemented("NGFRD::load(): TODO");
    }
    const Real volume() const override
    {
        const auto& edge = ps_->edge_lengths();
        return edge[0] * edge[1] * edge[2];
    }

    /**
     * return if the species is in this space or not.
     * @param sp a species
     * @return if the species is in this space
     */
    bool has_species(const Species& sp) const override
    {
        return ps_->has_species(sp);
    }

    std::vector<Species> list_species() const override
    {
        return ps_->list_species();
    }

    /**
     * get the number of molecules
     * @param sp a species
     * @return a number of molecules Integer
     */
    Integer num_molecules(const Species& sp) const override
    {
        return ps_->num_molecules(sp);
    }

    Integer num_molecules_exact(const Species& sp) const override
    {
        return ps_->num_molecules_exact(sp);
    }

    Real get_value(const Species& sp) const override
    {
        throw NotSupported("get_value(const Species&) is not supported"
                           " by this space class");
    }

    Real get_value_exact(const Species& sp) const override
    {
        throw NotSupported("get_value(const Species&) is not supported"
                           " by this space class");
    }

    /**
     * get the axes lengths of a cuboidal region.
     * @return edge lengths Real3
     */
    const Real3& edge_lengths() const override
    {
        return ps_->edge_lengths();
    }

    /**
     * get the number of particles.
     * @return a number of particles Integer
     */
    Integer num_particles() const override
    {
        return ps_->num_particles();
    }

    /**
     * get the number of particles.
     * @param sp a species
     * @return a number of particles Integer
     */
    Integer num_particles(const Species& sp) const override
    {
        return ps_->num_particles(sp);
    }

    Integer num_particles_exact(const Species& sp) const override
    {
        return ps_->num_particles_exact(sp);
    }

    /**
     * check if the particle exists.
     * @param pid an ID for the particle
     * @return if the particle exists or not bool
     */
    bool has_particle(const ParticleID& pid) const override
    {
        return ps_->has_particles(pid);
    }

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const override
    {
        return ps_->get_particles(pid);
    }

    /**
     * get all particles.
     * @return a list of particles
     */
    std::vector<std::pair<ParticleID, Particle>> list_particles() const
    {
        return ps_->list_particles();
    }

    /**
     * get particles.
     * @param sp a species
     * @return a list of particles
     */
    std::vector<std::pair<ParticleID, Particle>>
    list_particles(const Species& sp) const override
    {
        return ps_->list_particles(sp);
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const
    {
        return ps_->list_particles_exact(sp);
    }

private:

    template<typename Filter>
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_impl(const Real3& center, const Real radius) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> list;
        // TODO
        return list;
    }
    template<typename Filter>
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_impl(const std::pair<Real3, FaceID>& center,
            const Real radius) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> list;
        // TODO
        return list;
    }

private:

    std::shared_ptr<RandomNumberGenerator> rng_;
    std::weak_ptr<Model>                   model_;
    particle_id_generator_type             pidgen_;
    std::unique_ptr<particle_space_type>   ps_;
    std::shared_ptr<Polygon>               poly_;
    polygon_container_type                 poly_con_; // PID <-> FaceID

    // 2D gfrd specific cache
    Real estimated_possible_largest_particle_radius_;
    std::unordered_map<FaceID, std::array<ecell4::Segment, 6>> barriers_;

};

} // ngfrd
} // ecell4
#endif//ECELL4_NGFRD_NGFRD_WORLD_HPP
