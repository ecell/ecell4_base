#ifndef ECELL4_PARTICLE_SPACE_RTREE_IMPL_HPP
#define ECELL4_PARTICLE_SPACE_RTREE_IMPL_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/assert.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/index/predicates.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/version.hpp>

#if BOOST_VERSION >= 105400
#define  ECELL4_HAS_BOOST_STATIC_VECTOR 1
#include <boost/container/static_vector.hpp>
#endif//BOOST_VERSION

#if BOOST_VERSION >= 105800
#define  ECELL4_HAS_BOOST_SMALL_VECTOR 1
#include <boost/container/small_vector.hpp>
#endif//BOOST_VERSION

#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Integer3.hpp>

#ifdef WITH_HDF5
#include <ecell4/core/ParticleSpaceHDF5Writer.hpp>
#endif

#include <set>
#include <ciso646>

// enable boost to recognize ecell4::Real3 and ecell4::AABB {{{
namespace boost { namespace geometry { namespace traits {

// for ecell4::Real3 as boost.geometry's point concept

template<> struct tag<ecell4::Real3 > {typedef point_tag type;};
template<> struct dimension<ecell4::Real3> : boost::mpl::int_<3> {};
template<> struct coordinate_type<ecell4::Real3> {typedef ecell4::Real type;};
template<> struct coordinate_system<ecell4::Real3> {typedef boost::geometry::cs::cartesian type;};

template<> struct access<ecell4::Real3, 0>
{
    static inline double get(ecell4::Real3 const& p)             throw() {return p[0];}
    static inline void   set(ecell4::Real3& p, const ecell4::Real value) throw() {p[0] = value;}
};
template<> struct access<ecell4::Real3, 1>
{
    static inline double get(ecell4::Real3 const& p)             throw() {return p[1];}
    static inline void   set(ecell4::Real3& p, const ecell4::Real value) throw() {p[1] = value;}
};
template<> struct access<ecell4::Real3, 2>
{
    static inline double get(ecell4::Real3 const& p)             throw() {return p[2];}
    static inline void   set(ecell4::Real3& p, const ecell4::Real value) throw() {p[2] = value;}
};

// for ecell4::AABB as boost.geometry's box concept
// XXX NOTE: this concept requires AABB to be modifiable.

template<> struct tag<ecell4::AABB> {typedef box_tag type;};
template<> struct point_type<ecell4::AABB>
{
    typedef typename ecell4::Real3 type;
};

template <std::size_t D>
struct indexed_access<ecell4::AABB, min_corner, D>
{
    typedef typename coordinate_type<
        typename point_type<ecell4::AABB>::type>::type ct;

    static inline ct get(ecell4::AABB const& b) throw()
    {return geometry::get<D>(b.lower());}

    static inline void set(ecell4::AABB& b, ct const& value) throw()
    {geometry::set<D>(b.lower(), value);}
};
template <std::size_t D>
struct indexed_access<ecell4::AABB, max_corner, D>
{
    typedef typename coordinate_type<
        typename point_type<ecell4::AABB>::type>::type ct;

    static inline ct get(ecell4::AABB const& b) throw()
    {return geometry::get<D>(b.upper());}

    static inline void set(ecell4::AABB& b, ct const& value) throw()
    {geometry::set<D>(b.upper(), value);}
};

} /* traits */ } /* geometry */ }// boost }}}

namespace ecell4
{

class ParticleSpaceRTreeImpl
    : public ParticleSpace
{
public:
    typedef ParticleSpaceRTreeImpl self_type;

    // rtree: XXX the `rtree_algo_type` value should be tuned.
    typedef AABB box_type;
    typedef boost::tuple<box_type, ParticleID, Particle> rtree_value_type;
    typedef boost::geometry::index::quadratic<6, 2>      rtree_algo_type;
    typedef boost::geometry::index::rtree<rtree_value_type, rtree_algo_type>
            rtree_type;

    // container
    typedef ParticleSpace base_type;
    typedef ParticleSpace::particle_container_type particle_container_type;
    typedef utils::get_mapper_mf<ParticleID, particle_container_type::size_type>::type
            key_to_value_map_type;
    typedef std::pair<ParticleID, Particle>                  value_type;
    typedef typename particle_container_type::iterator       iterator;
    typedef typename particle_container_type::const_iterator const_iterator;

    // species support
    typedef std::set<ParticleID> particle_id_set;
    typedef utils::get_mapper_mf<Species::serial_type, particle_id_set>::type
            per_species_particle_id_set;

protected:

#ifdef ECELL4_HAS_BOOST_SMALL_VECTOR
    const static std::size_t typical_query_result_size = 10;
    const static std::size_t qsz = typical_query_result_size; // just an alias.

    typedef boost::container::small_vector<rtree_value_type, qsz>
            query_result_container_type;
#else
    typedef std::vector<rtree_value_type> query_result_container_type;
#endif//ECELL4_HAS_BOOST_SMALL_VECTOR


#ifdef ECELL4_HAS_BOOST_STATIC_VECTOR
    typedef boost::container::static_vector<box_type, 8>
            query_boxes_container_type;
#else
    typedef std::vector<box_type> query_boxes_container_type;
#endif//ECELL4_HAS_BOOST_STATIC_VECTOR


    struct particle_id_excluder
    {
        ParticleID pid;

        explicit particle_id_excluder(const ParticleID& ignore): pid(ignore){}

        inline bool operator()(const rtree_value_type& x) const
        {
            return boost::get<1>(x) != pid;
        }
    };

    struct particle_id2_excluder
    {
        ParticleID pid1, pid2;

        particle_id2_excluder(const ParticleID& ignore1, const ParticleID& ignore2)
            : pid1(ignore1), pid2(ignore2)
        {}

        inline bool operator()(const rtree_value_type& x) const
        {
            return (boost::get<1>(x) != pid1) && (boost::get<1>(x) != pid2);
        }
    };

public:

    explicit ParticleSpaceRTreeImpl(const Real3& edge_lengths)
        : base_type(), max_radius_(0.0), edge_lengths_(edge_lengths)
    {}

    void reset(const Real3& edge_lengths);

    // inherit from Space

    virtual Integer num_species() const
    {
        return particle_pool_.size();
    }

    virtual bool has_species(const Species& sp) const
    {
        return (particle_pool_.find(sp.serial()) != particle_pool_.end());
    }

    virtual std::vector<Species> list_species() const;

    // inherit from ParticleSpace

    const Real3& edge_lengths() const
    {
        return edge_lengths_;
    }

    const particle_container_type& particles() const
    {
        return particles_;
    }

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {
        const particle_container_type::const_iterator found = this->find(pid);
        if(found == this->particles_.end())
        {
            throw NotFound("ParticleSpaceRTree::get_particle: particle not found.");
        }
        return *found;
    }

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        if(this->has_particle(pid))
        {
            this->update(pid, p);
            return false;
        }
        else
        {
            this->insert(pid, p);
            return true;
        }
    }
    bool has_particle(const ParticleID& pid) const
    {
        return (this->find(pid) != this->particles_.end());
    }
    void remove_particle(const ParticleID& pid)
    {
        this->remove(pid);
    }

    Integer num_particles() const
    {
        return particles_.size();
    }
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    std::vector<std::pair<ParticleID, Particle> >
    list_particles() const
    {
        return particles_;
    }
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const;

    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        save_particle_space(*this, root);
    }

    void load_hdf5(const H5::Group& root)
    {
        load_particle_space(root, this);
    }
#endif

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius,
            const ParticleID& ignore) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const;

    // this requires boost::geometry::index::* as query.
    template<typename Query>
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_satisfies(Query query) const
    {
        query_result_container_type tmp;
        rtree_.query(query, std::back_inserter(tmp));

        std::vector<std::pair<ParticleID, Particle> > retval;
        retval.reserve(tmp.size());
        for(typename query_result_container_type::const_iterator
                i(tmp.begin()), e(tmp.end()); i!=e; ++i)
        {
            retval.push_back(
                    std::make_pair(boost::get<1>(*i), boost::get<2>(*i)));
        }
        return retval;
    }

protected:

    // XXX: this assumes ecell4::Particle has spherical shape.
    static inline box_type make_box(const Particle& p)
    {
        const Real r = p.radius();
        Real3 up(p.position()), lw(p.position());
        up[0] += r; up[1] += r; up[2] += r;
        lw[0] -= r; lw[1] -= r; lw[2] -= r;
        return box_type(lw, up);
    }
    static inline box_type make_box(const Real3& center, const Real radius)
    {
        Real3 up(center), lw(center);
        up[0] += radius; up[1] += radius; up[2] += radius;
        lw[0] -= radius; lw[1] -= radius; lw[2] -= radius;
        return box_type(lw, up);
    }

    template<std::size_t N>
    inline typename boost::enable_if_c<(N<3), void>::type
    split_box_by_boundary(query_boxes_container_type& boxes) const
    {
        // assuming the boundary is [0, edge_lengths_)
        //          and the query box is inside of the boundary...
        const std::size_t sz = boxes.size();
        for(std::size_t i=0; i<sz; ++i)
        {
            if(this->edge_lengths_[N] < boxes[i].upper()[N])
            {
                box_type bx(boxes[i]);
                bx.lower()[N] = 0.0;
                bx.upper()[N] = boxes[i].upper()[N] - this->edge_lengths_[N];
                boxes[i].upper()[N] = this->edge_lengths_[N];
                boxes.push_back(bx);// XXX do not use iterator
            }
            else if(boxes[i].lower()[N] < 0)
            {
                box_type bx(boxes[i]);
                bx.lower()[N] = boxes[i].lower()[N] + this->edge_lengths_[N];
                bx.upper()[N] = this->edge_lengths_[N];
                boxes[i].lower()[N] = 0.0;
                boxes.push_back(bx);// XXX do not use iterator
            }
        }
        return;
    }

    static inline rtree_value_type
    make_rtree_value(const ParticleID& pid, const Particle& p)
    {
        return boost::make_tuple(make_box(p), pid, p);
    }

    inline particle_container_type::iterator
    find(const ParticleID& k)
    {
        key_to_value_map_type::const_iterator p(idx_map_.find(k));
        return (idx_map_.end() == p) ? particles_.end() :
                                       particles_.begin() + p->second;
    }

    inline particle_container_type::const_iterator
    find(const ParticleID& k) const
    {
        key_to_value_map_type::const_iterator p(idx_map_.find(k));
        return (idx_map_.end() == p) ? particles_.end() :
                                       particles_.begin() + p->second;
    }

    void insert(const ParticleID& pid, const Particle& p)
    {
        if(this->has_particle(pid))
        {
            throw std::invalid_argument("ParticleSpaceRTree::insert: already has");
        }
        rtree_.insert(make_rtree_value(pid, p));
        particle_pool_[p.species_serial()].insert(pid);
        if(max_radius_ < p.radius())
        {
            max_radius_ = p.radius();
        }

        const std::size_t idx = particles_.size();
        particles_.push_back(std::make_pair(pid, p));
        idx_map_[pid] = idx;
    }

    void remove(const ParticleID& pid)
    {
        if(not this->has_particle(pid))
        {
            throw NotFound("ParticleSpaceRTree::remove: particle not found");
        }
        const std::size_t idx = idx_map_[pid];
        idx_map_.erase(pid);

        const Particle& p = this->particles_[idx].second;
        const std::size_t result = rtree_.remove(make_rtree_value(pid, p));
        BOOST_ASSERT(result == 1);

        particle_pool_[p.species_serial()].erase(pid);

        idx_map_[particles_.back().first] = idx;
        particles_[idx] = particles_.back();
        particles_.pop_back();

        return;
    }

    void update(const ParticleID& pid, const Particle& p)
    {
        if(not this->has_particle(pid))
        {
            throw NotFound("ParticleSpaceRTree::update: particle not found");
        }
        const std::size_t idx = idx_map_[pid];
        const std::size_t result = rtree_.remove(make_rtree_value(pid, particles_[idx].second));
        BOOST_ASSERT(result == 1);
        rtree_.insert(make_rtree_value(pid, p));

        if(particles_[idx].second.species() != p.species())
        {
            particle_pool_[particles_[idx].second.species_serial()].erase(pid);
            particle_pool_[p.species_serial()].insert(pid);

            if(max_radius_ < p.radius())
            {
                max_radius_ = p.radius();
            }
        }
        particles_[idx].second = p;
        return;
    }

protected:

    Real  max_radius_; // XXX
    Real3 edge_lengths_;

    rtree_type                  rtree_;
    particle_container_type     particles_;
    key_to_value_map_type       idx_map_;
    per_species_particle_id_set particle_pool_;
};

}; // ecell4

#endif /* ECELL4_PARTICLE_SPACE_CELL_LIST_IMPL_HPP */
