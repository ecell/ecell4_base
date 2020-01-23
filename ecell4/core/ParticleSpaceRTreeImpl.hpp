#ifndef ECELL4_PARTICLE_SPACE_RTREE_IMPL_HPP
#define ECELL4_PARTICLE_SPACE_RTREE_IMPL_HPP

#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Integer3.hpp>
#include <ecell4/core/PeriodicRTree.hpp>

#ifdef WITH_HDF5
#include <ecell4/core/ParticleSpaceHDF5Writer.hpp>
#endif

#include <set>

namespace ecell4
{

class ParticleSpaceRTreeImpl
    : public ParticleSpace
{
public:

    struct ParticleAABBGetter
    {
        AABB operator()(const Particle& p, const Real margin) const noexcept
        {
            const Real3 radius(p.radius() + p.D() * margin,
                               p.radius() + p.D() * margin,
                               p.radius() + p.D() * margin);
            return AABB(p.position() - radius, p.position() + radius);
        }
    };

    using base_type  = ParticleSpace;
    using rtree_type = PeriodicRTree<ParticleID, Particle, ParticleAABBGetter>;
    using box_type                = typename rtree_type::box_type;
    using value_type              = typename rtree_type::value_type;
    using key_to_value_map_type   = typename rtree_type::key_to_value_map_type;
    using particle_container_type = typename rtree_type::container_type;
    using iterator                = typename rtree_type::iterator;
    using const_iterator          = typename rtree_type::const_iterator;

    static_assert(std::is_same<value_type,
            std::pair<ParticleID, Particle>>::value, "");
    static_assert(std::is_same<particle_container_type,
            ParticleSpace::particle_container_type>::value, "");

    // species support
    using particle_id_set             = std::set<ParticleID>;
    using per_species_particle_id_set =
        utils::get_mapper_mf<Species::serial_type, particle_id_set>::type;

public:

    // the default value of margin should be tuned later.
    explicit ParticleSpaceRTreeImpl(const Real3& edge_lengths,
                                    const Real margin = 0.1)
        : base_type(), rtree_(edge_lengths, margin)
    {}

    void reset(const Real3& edge_lengths);

    // inherit from Space

    virtual Integer num_species() const
    {
        return particle_pool_.size();
    }

    virtual bool has_species(const Species& sp) const
    {
        return (particle_pool_.count(sp.serial()) != 0);
    }

    // inherit from ParticleSpace

    const Real3& edge_lengths() const override
    {
        return rtree_.edge_lengths();
    }

    const particle_container_type& particles() const override
    {
        // Since PeriodicRTree manages particles as their index, the container
        // can occasionally contains "already-erased" particles marked as
        // "overwritable", like a colony (P0447R1). Because of this, the raw
        // reference may contain an invalid particle that is already removed
        // from the space. It may invalidates subsequent operation.
        //     To avoid confusion, this method throws an exception. This will
        // never be implemented, so it throws a `NotSupported`, not a
        // `NotImplemented`. Since the default implementation of `list_species()`
        // uses this function, it is overriden later.
        throw NotSupported("ParticleSpaceRTreeImpl does not support "
                           "`particle_container_type const& particles()`");
    }

    // ParticleSpace has the default list_species implementation.
    // But it uses particles() member method that cannot be implemented with
    // PeriodicRTree. So it overwrites the default implementation.
    std::vector<Species> list_species() const override;

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {
        if(!rtree_.has(pid))
        {
            throw NotFound("ParticleSpaceRTreeImpl::get_particle(): "
                           "particle not found.");
        }
        return rtree_.get(pid);
    }

    bool update_particle(const ParticleID& pid, const Particle& newp)
    {
        if(rtree_.has(pid))
        {
            const auto& oldp = rtree_.get(pid).second;
            if(oldp.species() != newp.species())
            {
                particle_pool_[oldp.species_serial()].erase (pid);
                particle_pool_[newp.species_serial()].insert(pid);
            }
        }
        particle_pool_[newp.species_serial()].insert(pid);
        return rtree_.update(pid, newp);
    }
    bool has_particle(const ParticleID& pid) const
    {
        return rtree_.has(pid);
    }
    void remove_particle(const ParticleID& pid)
    {
        const auto& p = rtree_.get(pid).second;
        particle_pool_[p.species_serial()].erase(pid);
        rtree_.erase(pid, p);
        return;
    }

    Integer num_particles() const
    {
        return rtree_.size();
    }
    Integer num_particles      (const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;
    Integer num_molecules      (const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    std::vector<std::pair<ParticleID, Particle> >
    list_particles() const
    {
        return rtree_.list_objects();
    }
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const;

    virtual void save(const std::string& filename) const
    {
        throw NotImplemented(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        throw NotImplemented(
            "save(const std::string) is not supported by this space class");
    }

    void load_hdf5(const H5::Group& root)
    {
        throw NotImplemented(
            "save(const std::string) is not supported by this space class");
    }
#endif

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& pos, const Real& radius) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& pos, const Real& radius,
            const ParticleID& ignore) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const;

protected:

    template<typename Query, typename OutputIterator>
    void query_impl(Query&& q, OutputIterator out) const
    {
        rtree_.query(std::forward<Query>(q), out);
        return ;
    }

    // ------------------------------------------------------------------------
    // query objects

    template<typename Filter>
    struct IntersectionQuery
    {
        Real3  center;
        Real   radius;
        Filter ignores;

        IntersectionQuery(const Real3& c, const Real r, Filter f) noexcept
            : center(c), radius(r), ignores(std::move(f))
        {}
        IntersectionQuery(const IntersectionQuery&) = default;
        IntersectionQuery(IntersectionQuery&&)      = default;
        IntersectionQuery& operator=(const IntersectionQuery&) = default;
        IntersectionQuery& operator=(IntersectionQuery&&)      = default;
        ~IntersectionQuery() = default;


        // If it does not matches, return boost::none.
        // If it matches, return pairof(pidp, distance).
        boost::optional<std::pair<value_type, Real>>
        operator()(const value_type& pidp, const Real3& edges) const noexcept
        {
            if(ignores(pidp)){return boost::none;}

            const auto rr  = this->radius + pidp.second.radius();
            const auto rhs = this->periodic_transpose(pidp.second.position(),
                                                      center, edges);
            const auto dist_sq = length_sq(rhs - center);
            if(rr * rr < dist_sq)
            {
                return boost::none;
            }
            return std::make_pair(pidp, std::sqrt(dist_sq));
        }

        bool operator()(const AABB& box, const Real3& edges) const noexcept
        {
            return this->distance_sq(box, this->center, edges) <
                   this->radius * this->radius;
        }

        // -------------------------------------------------------------------
        // geometry stuffs

        // AABB-sphere intersection query under the PBC
        Real distance_sq(const AABB& box, Real3 pos, const Real3& edge_lengths) const noexcept
        {
            pos = periodic_transpose(pos, (box.upper() + box.lower()) * 0.5, edge_lengths);

            Real dist_sq = 0;
            for(std::size_t i=0; i<3; ++i)
            {
                const auto v = pos[i];
                if(v < box.lower()[i])
                {
                    dist_sq += (v - box.lower()[i]) * (v - box.lower()[i]);
                }
                else if(box.upper()[i] < v)
                {
                    dist_sq += (v - box.upper()[i]) * (v - box.upper()[i]);
                }
            }
            return dist_sq;
        }

        // transpose a position based on the periodic boundary condition.
        Real3 periodic_transpose(
            const Real3& pos1, const Real3& pos2, const Real3& edges) const
        {
            Real3 retval(pos1);
            for(std::size_t dim(0); dim < 3; ++dim)
            {
                const Real edge_length(edges[dim]);
                const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

                if (half < diff)
                {
                    retval[dim] += edge_length;
                }
                else if (diff < -half)
                {
                    retval[dim] -= edge_length;
                }
            }
            return retval;
        }
    };

    template<typename Filter>
    static IntersectionQuery<Filter> make_intersection_query(
            const Real3& c, const Real r, Filter&& f)
    {
        return IntersectionQuery<Filter>(c, r, std::forward<Filter>(f));
    }

protected:

    rtree_type                  rtree_;
    per_species_particle_id_set particle_pool_;
};

}; // ecell4

#endif /* ECELL4_PARTICLE_SPACE_CELL_LIST_IMPL_HPP */
