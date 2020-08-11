#ifndef ECELL4_NGFRD_NGFRD_WORLD_HPP
#define ECELL4_NGFRD_NGFRD_WORLD_HPP

#include <memory>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Segment.hpp>
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

//
// Since most of the member methods are redirected to ParticleSpace, we store
// both 3D and 2D particles to the same ParticleSpace, as spherical particles.
// It causes redundant checking in some cases because 2D particles which have
// circular shape are considered sphere in ParticleSpace class. We need to check
// if the particle is on a face and the target collides with the face before
// concluding particles collide with each other.
//
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
               const Integer3& matrix_sizes = Integer3(3, 3, 3),
               const Real      margin       = 0.1)
        : ps_(new particle_space_type(edge_lengths, margin)),
          polygon_(std::make_shared<Polygon>(edge_lengths, matrix_sizes))
    {
        rng_ = std::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();

        this->prepare_restrictions();
    }

    NGFRDWorld(const Real3& edge_lengths, const Integer3& matrix_sizes,
               const Real margin, // bounding box margin
               std::shared_ptr<RandomNumberGenerator> rng,
               std::shared_ptr<Polygon>               poly)
        : ps_(new particle_space_type(edge_lengths, margin)),
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
    new_particle(const Particle& p)
    {
        ParticleID pid(this->pidgen_());

        if( ! has_overlapping_faces(p.position(), p.radius()).empty() &&
           this->list_particles_within_radius_3D(p.position(), p.radius()))
        {
            // XXX: do NOT call this->update_particle to avoid redundant check
            this->ps_->update_particle(pid, p);
            return std::make_pair(std::make_pair(pid, p), true);
        }
        else
        {
            return std::make_pair(std::make_pair(pid, p), false);
        }
    }

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p, const FaceID& fid)
    {
        ParticleID pid(this->pidgen_());

        if(this->list_particles_within_radius_2D(
                    std::make_pair(p.position(), fid), p.radius()).empty())
        {
            this->ps_->update_particle(pid, p);
            this->poly_con_.update(pid, fid);
            return std::make_pair(std::make_pair(pid, p), true);
        }
        else
        {
            return std::make_pair(std::make_pair(pid, p), false);
        }
    }

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        if( ! has_overlapping_faces(p.position(), p.radius()).empty() &&
           this->list_particles_within_radius_3D(p.position(), p.radius()))
        {
            // if `pid` already exists and was a 2D particle, we need to reset
            // the relationship.
            if(const auto fid = this->poly_con_.on_which_face(pid))
            {
                this->poly_con_.remove(pid, *fid);
            }
            return this->ps_->update_particle(pid, p);
        }
        else
        {
            // overlap found. the update is rejected. no change.
            return true;
        }
    }
    bool update_particle(const ParticleID& pid, const Particle& p,
                         const FaceID& fid)
    {
        if(this->list_particles_within_radius_2D(
                    std::make_pair(p.position(), fid), p.radius()).empty())
        {
            this->ps_->update_particle(pid, p);
            this->poly_con_.update(pid, fid);
            return std::make_pair(std::make_pair(pid, p), true);
        }
        else
        {
            return std::make_pair(std::make_pair(pid, p), false);
        }
    }

    void remove_particle(const ParticleID& pid)
    {
        if(const auto fid = poly_con_.on_which_face(pid))
        {
            poly_con_.remove(pid, *fid);
        }
        ps_->remove_particle(pid);
        return;
    }
    bool has_particle(const ParticleID& pid)
    {
        return ps_->has_particle(pid);
    }

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
    // for speedup

    bool has_overlapping_particles_3D(const Real3& center, const Real radius) const
    {
        return !(this->list_particles_within_radius_3D(center, radius).empty());
    }
    bool has_overlapping_particles_3D(const Real3& center, const Real radius,
            const ParticleID& ignore) const
    {
        return !(this->list_particles_within_radius_3D(center, radius, ignore).empty());
    }
    bool has_overlapping_particles_3D(const Real3& center, const Real radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return !(this->list_particles_within_radius_3D(
                    center, radius, ignore1, ignore2).empty());
    }

    bool has_overlapping_particles_2D(
            const std::pair<Real3, FaceID>& center, const Real radius) const
    {
        return !(this->list_particles_within_radius_2D(center, radius).empty());
    }
    bool has_overlapping_particles_2D(
            const std::pair<Real3, FaceID>& center, const Real radius,
            const ParticleID& ignore) const
    {
        return !(this->list_particles_within_radius_2D(
                    center, radius, ignore).empty());
    }
    bool has_overlapping_particles_2D(
            const std::pair<Real3, FaceID>& center, const Real radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return !(this->list_particles_within_radius_2D(
                    center, radius, ignore1, ignore2).empty());
    }

    bool has_overlapping_faces(const Real3& center, const Real radius) const
    {
        return polygon_->has_overlapping_faces(center, radius);
    }
    bool has_overlapping_faces(const Real3& center, const Real radius,
                               const FaceID& ignore) const
    {
        return polygon_->has_overlapping_faces(center, radius, ignore);
    }

    // ------------------------------------------------------------------------
    // list_particles_within_radius_3D
    // Note that it returns only 3D particles. No 2D-3D checking.

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_3D(const Real3& center, const Real radius) const
    {
        return this->list_particles_within_radius_3D_impl(center, radius);
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_3D(const Real3& center, const Real radius,
            const ParticleID& ignore) const
    {
        return this->list_particles_within_radius_3D_impl(
                center, radius, ignore);
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_3D(const Real3& center, const Real radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return this->list_particles_within_radius_3D_impl(
                center, radius, ignore1, ignore2);
    }

    // ------------------------------------------------------------------------
    // list_particles_within_radius_2D
    // Note that it returns only 2D particles. No 2D-3D checking.

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_2D(
            const std::pair<Real3, FaceID>& center, const Real radius) const
    {
        return this->list_particles_within_radius_2D_impl(center, radius,
                [](const ParticleID& pid) noexcept -> bool {return false;});
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_2D(
            const std::pair<Real3, FaceID>& center, const Real radius,
            const ParticleID& ignore) const
    {
        return this->list_particles_within_radius_2D_impl(center, radius,
                [](const ParticleID& pid) noexcept -> bool {
                    return pid == ignore;
                });
    }
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_2D(
            const std::pair<Real3, FaceID>& center, const Real radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return this->list_particles_within_radius_2D_impl(center, radius,
                [](const ParticleID& pid) noexcept -> bool {
                    return pid == ignore1 || pid == ignore2;
                });
    }

    Polygon&       polygon()       noexcept {return *polygon_;}
    Polygon const& polygon() const noexcept {return *polygon_;}

    PeriodicBoundary boundary() const noexcept
    {
        return PeriodicBoundary(this->edge_lengths());
    }

    molecule_info_type get_molecule_info(const Species& sp) const
    {
        Real radius(0.0), D(0.0);

        if (sp.has_attribute("radius") && sp.has_attribute("D"))
        {
            radius = sp.get_attribute_as<Real>("radius");
            D = sp.get_attribute_as<Real>("D");
        }
        else if (std::shared_ptr<Model> bound_model = lock_model())
        {
            Species newsp(bound_model->apply_species_attributes(sp));
            if (newsp.has_attribute("radius")
                && newsp.has_attribute("D"))
            {
                radius = newsp.get_attribute_as<Real>("radius");
                D = newsp.get_attribute_as<Real>("D");
            }
        }

        if (radius <= 0.0)
        {
            throw_exception<IllegalArgument>("A particle with invalid size [",
                    radius, "] was given.");
        }

        return MoleculeInfo{radius, D};
    }

    void bind_to(std::shared_ptr<Model> model)
    {
        if (std::shared_ptr<Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to NGFRDWorld"
                    << std::endl;
            }
        }
        model_ = model;
    }

    std::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
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

    template<typename ... Ts>
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_impl(const Real3& center, const Real radius,
                                      Ts&& ... ignores) const
    {
        const auto faces = polygon_->list_faces_within_radius(center, radius);
        auto list = this->ps_->list_particles_within_radius(
                center, radius, std::forward<Ts>(ignores)...);

        // FIXME
        // In case of a particle straddles faces, the calculation becomes
        // incorrect (overestimate the number of particles). We need to remove
        // them to make it perfect.
        //
        // To correctly find the overlap, we first need to check if a particle
        // overlaps with a face. If
        // - the range overlaps with a face and
        // - a particle overlaps with the face and
        // - the distance between the center of the particle and the center of
        //   the region is less than the radius and the range radius on the face
        // then they overlaps each other.
        //
        // If particle does not straddle faces, it may not overlap with the
        // region even if the distance between region and the particle is less
        // than the radius.

        // 3D particle collides with 2D particle only if it is within radius
        // AND 3D particle collides with the face on which the 2D particle is.

        const auto removed = std::remove_if(list.begin(), list.end(), [this, &faces]
            (const std::pair<std::pair<ParticleID, Particle>, Real>& elem) -> bool
            {
                const auto pid = elem.first.first;
                if(const auto fidopt = poly_con_.on_which_face(elem.first.first))
                {
                    // 2D particle! If the face is not found in the overlapping
                    // face, remove it.
                    const auto fid = *fidopt;

                    return std::find_if(faces.begin(), faces.end(), [&fid]
                        (const std::pair<std::pair<FaceID, Triangle>, Real>& f) {
                            return fid == f.first.first;
                        }) == faces.end();
                }
                // it is not on a face. it is 3D. do not remove this.
                return false;
            });

        // erase all the 2D particles that do not collide actually.
        list.erase(removed, list.end());
        return list;
    }


    template<typename ... Ts>
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_3D_impl(const Real3& center, const Real radius,
                                         Ts&& ... ignores) const
    {
        auto list = this->ps_->list_particles_within_radius(
                center, radius, std::forward<Ts>(ignores)...);

        // It does not check collision with 2D particles.
        // If the list contains 2D particle, remove it from the result.
        const auto removed = std::remove_if(list.begin(), list.end(),
            [this](const std::pair<std::pair<ParticleID, Particle>, Real>& elem)
                noexcept -> bool {
                return poly_con_.on_face(elem.first.first);
            });

        // erase all the 2D particles that do not collide actually.
        list.erase(removed, list.end());
        return list;
    }

    template<typename Filter>
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius_2D_impl(
            const std::pair<Real3, FaceID>& center, const Real radius,
            Filter filter) const
    {
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> list;

        const auto& pos = center.first;
        const auto& fid = center.second;

        // check particles on the same face
        for(const auto& pid : this->poly_con_.objects_on(fid))
        {
            auto pp = ps_->get_particle(pid);
            // it is okay to use 3D distance because both are on the same face
            const Real dist = length(pos - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(std::move(pp), dist));
            }
        }

        // check particles on the neighborling faces
        for(const auto& fid : polygon_->neighbor_faces_of(fid))
        {
            for(const auto& pid : this->poly_con_.objects_on(fid))
            {
                auto pp = ps_->get_particle(pid);
                const Real dist = ecell4::polygon::distance(*polygon_,
                    pos, std::make_pair(pp.second.position(), fid)
                    ) - pp.second.radius();

                if(dist < radius)
                {
                    retval.push_back(std::make_pair(std::move(pp), dist));
                }
            }
        }

        // no 3D check (Since particle would never overlaps with any face).
        return list;
    }

  private:

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
                const auto pvi1 = this->boundary().periodic_transpose(
                        polygon_->position_at(vi1), pvi0);
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

    std::unique_ptr<particle_space_type>   ps_;
    std::shared_ptr<RandomNumberGenerator> rng_;
    std::weak_ptr<Model>                   model_;
    SerialIDGenerator<ParticleID>          pidgen_;
    std::shared_ptr<Polygon>               polygon_;
    polygon_container_type                 poly_con_; // PID <-> FaceID

    // 2D gfrd specific cache
    Real estimated_possible_largest_particle_radius_;
    std::unordered_map<FaceID, std::array<ecell4::Segment, 6>> barriers_;

};

} // ngfrd
} // ecell4
#endif//ECELL4_NGFRD_NGFRD_WORLD_HPP
