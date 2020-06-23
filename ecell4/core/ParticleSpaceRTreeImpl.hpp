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
        std::unordered_map<Species::serial_type, particle_id_set>;

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

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const override
    {
        if(!rtree_.has(pid))
        {
            std::ostringstream oss;
            oss << "ParticleSpaceRTreeImpl::get_particle: No such particle ("
                << pid << ").";
            throw NotFound(oss.str());
        }
        return rtree_.get(pid);
    }

    // returns true if it adds a new particle
    bool update_particle(const ParticleID& pid, const Particle& newp) override
    {
        if(rtree_.has(pid))
        {
            const auto& oldp = rtree_.get(pid).second;
            if(oldp.species() != newp.species())
            {
                particle_pool_[oldp.species_serial()].erase (pid);
                particle_pool_[newp.species_serial()].insert(pid);
            }
            // if species does not change, then we don't need to do anything.
        }
        else
        {
            // if `newp` is completely new, we need to insert it to the pool.
            particle_pool_[newp.species_serial()].insert(pid);
        }
        const auto retval = rtree_.update(pid, newp);
        assert(rtree_.diagnosis());
        assert(this->diagnosis());
        return retval;
    }
    bool has_particle(const ParticleID& pid) const override
    {
        return rtree_.has(pid);
    }
    void remove_particle(const ParticleID& pid) override
    {
        if(!rtree_.has(pid))
        {
            std::ostringstream oss;
            oss << "ParticleSpaceRTreeImpl::remove_particle: No such particle ("
                << pid << ").";
            throw NotFound(oss.str());
        }
        const auto& p = rtree_.get(pid).second;
        particle_pool_[p.species_serial()].erase(pid);
        rtree_.erase(pid, p);
        return;
    }

    Integer num_particles() const override
    {
        return rtree_.size();
    }
    Integer num_particles      (const Species& sp) const override;
    Integer num_particles_exact(const Species& sp) const override;
    Integer num_molecules      (const Species& sp) const override;
    Integer num_molecules_exact(const Species& sp) const override;

    std::vector<std::pair<ParticleID, Particle> >
    list_particles() const override
    {
        return rtree_.list_objects();
    }
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const override;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const override;

    virtual void save(const std::string& filename) const
    {
        throw NotImplemented(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const override
    {
        save_particle_space(*this, root);
    }

    void load_hdf5(const H5::Group& root) override
    {
        load_particle_space(root, this);
    }
#endif

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& pos, const Real& radius) const override;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& pos, const Real& radius,
            const ParticleID& ignore) const override;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
    list_particles_within_radius(const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const override;

    bool diagnosis() const
    {
        bool is_ok = true;
        const auto ps = this->list_particles();
        for(std::size_t i=0; i+1<ps.size(); ++i)
        {
            const auto& pi = ps.at(i).second;
            for(std::size_t j=i+1; j<ps.size(); ++j)
            {
                const auto& pj = ps.at(j).second;
                const auto d2 = this->distance_sq(pi.position(), pj.position());
                const auto dist = std::sqrt(d2) - pi.radius() - pj.radius();

                if(dist < 0.0)
                {
                    std::cerr << "particle " << ps.at(i).first;
                    std::cerr << " and " << ps.at(j).first;
                    std::cerr << " collide with each other" << std::endl;
                    is_ok = false;
                }
            }
        }
        return is_ok;
    }

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
        operator()(const value_type& pidp, const PeriodicBoundary& pbc) const noexcept
        {
            if(ignores(pidp)){return boost::none;}

            // use the same algorithm as the ParticleSpaceVectorImpl.
            const auto rhs = pbc.periodic_transpose(pidp.second.position(),
                                                    this->center);
            const auto dist = length(this->center - rhs) - pidp.second.radius();
            if(dist <= this->radius)
            {
                return std::make_pair(pidp, dist);
            }
            return boost::none;
        }

        bool operator()(const AABB& box, const PeriodicBoundary& pbc) const noexcept
        {
            return this->distance_sq(box, this->center, pbc) <=
                   this->radius * this->radius;
        }

        // -------------------------------------------------------------------
        // AABB-sphere distance under the PBC
        Real distance_sq(const AABB& box, Real3 pos, const PeriodicBoundary& pbc) const noexcept
        {
            pos = pbc.periodic_transpose(pos, (box.upper() + box.lower()) * 0.5);

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
