#ifndef MULTI_HPP
#define MULTI_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>

#include <ecell4/core/functions.hpp>
using ecell4::pow_2;
#include <ecell4/core/comparators.hpp>
#ifdef WIN32_MSC
#include <boost/container/map.hpp>
#endif
namespace win32_workaround
{
template <typename Key, typename T>
struct get_mapper_mf
{
#ifdef WIN32_MSC
    typedef boost::container::map<Key, T> type;
#else
    typedef std::map<Key, T> type;
#endif
};
}

#include "exceptions.hpp"
#include "Domain.hpp"
#include "ParticleContainer.hpp"
#include "ParticleContainerBase.hpp"
// #include "Sphere.hpp"
// #include "BDSimulator.hpp"
#include "BDPropagator.hpp"
#include "Logger.hpp"
#include "VolumeClearer.hpp"
#include "utils/array_helper.hpp"
#include "utils/range.hpp"

#include <greens_functions/PairGreensFunction.hpp>
// using namespace greens_functions;

template<typename Ttraits_>
class MultiParticleContainer
    : public Ttraits_::world_type::particle_container_type
{
public:

    typedef typename Ttraits_::world_type::particle_container_type base_type;
    typedef typename Ttraits_::world_type world_type;
    typedef typename world_type::traits_type traits_type;
    // typedef ParticleContainerUtils<traits_type> utils;

    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename traits_type::molecule_info_type molecule_info_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::particle_id_pair_generator
        particle_id_pair_generator;
    typedef typename traits_type::particle_id_pair_and_distance
        particle_id_pair_and_distance;
    typedef typename traits_type::particle_id_pair_and_distance_list
        particle_id_pair_and_distance_list;

    typedef typename base_type::transaction_type transaction_type;

    //XXX: typedef std::map<particle_id_type, particle_type> particle_map;
    typedef typename win32_workaround::get_mapper_mf<particle_id_type, particle_type>::type particle_map;
    typedef sized_iterator_range<typename particle_map::const_iterator> particle_id_pair_range;

    typedef typename world_type::particle_container_type::time_type time_type;

    virtual ~MultiParticleContainer() {}

    virtual ecell4::Integer num_particles() const
    {
        return particles_.size();
    }

    // virtual size_type num_particles() const
    // {
    //     return particles_.size();
    // }

    virtual const position_type& edge_lengths() const
    {
        return world_.edge_lengths();
    }

    // virtual molecule_info_type const& find_molecule_info(species_id_type const& id) const
    // {
    //     return world_.find_molecule_info(id);
    // }

    // virtual molecule_info_type const& get_molecule_info(species_id_type const& id)
    // {
    //     return world_.get_molecule_info(id);
    // }

    virtual molecule_info_type get_molecule_info(species_id_type const& id) const
    {
        return world_.get_molecule_info(id);
    }

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return world_.get_structure(id);
    }

    virtual std::pair<particle_id_pair, bool> new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        std::pair<particle_id_pair, bool> const retval(world_.new_particle(sid, pos));
        particles_.insert(retval.first);
        return retval;
    }

    virtual bool update_particle(const particle_id_type& pid, const particle_type& p)
    {
        world_.update_particle(pid, p);
        typename particle_map::iterator const i(particles_.find(pid));
        if (i != particles_.end())
        {
            (*i).second = p;
            return false;
        }
        else
        {
            particles_.insert(i, std::make_pair(pid, p));
            return true;
        }
    }

    virtual void remove_particle(particle_id_type const& id)
    {
        world_.remove_particle(id);
        particles_.erase(id);
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        typename particle_map::const_iterator i(particles_.find(id));
        if (particles_.end() == i)
        {
            throw not_found(std::string("No such particle: id=")
                    + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return particles_.end() != particles_.find(id);
    }

    virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s) const
    {
        return check_overlap(s, array_gen<particle_id_type>());
    }

    virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return check_overlap(s, array_gen(ignore));
    }

    virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return check_overlap(s, array_gen(ignore1, ignore2));
    }

    template<typename Tsph_, typename Tset_>
    particle_id_pair_and_distance_list check_overlap(Tsph_ const& s, Tset_ const& ignore) const
    {
        particle_id_pair_and_distance_list retval;
        for (typename particle_map::const_iterator i(particles_.begin()),
                                                   e(particles_.end());
             i != e; ++i)
        {
            length_type const dist(world_.distance(shape((*i).second), s.position()));
            if (dist < s.radius() && !collection_contains(ignore, (*i).first))
            {
                retval.push_back(std::make_pair(*i, dist));
            }
        }
        std::sort(retval.begin(), retval.end(),
            ecell4::utils::pair_second_element_comparator<particle_id_pair, length_type>());
        return retval;
    }

    virtual transaction_type* create_transaction()
    {
        return new TransactionImpl<MultiParticleContainer>(*this);
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return world_.distance(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return world_.apply_boundary(v);
    }

    // virtual length_type apply_boundary(length_type const& v) const
    // {
    //     return world_.apply_boundary(v);
    // }

    virtual position_type periodic_transpose(position_type const& p0, position_type const& p1) const
    {
        return world_.periodic_transpose(p0, p1);
    }

    // virtual length_type periodic_transpose(length_type const& p0, length_type const& p1) const
    // {
    //     return world_.periodic_transpose(p0, p1);
    // }

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range(particles_.begin(), particles_.end(),
                                      particles_.size());
    }

    MultiParticleContainer(world_type& world): world_(world) {}

    /** ecell4::Space
     */
    // virtual const time_type& t() const
    virtual const time_type t() const
    {
        return world_.t();
    }

    virtual void set_t(const time_type& t)
    {
        world_.set_t(t);
    }

private:
    world_type& world_;
    particle_map particles_;
};

template<typename Tsim_>
class Multi: public Domain<typename Tsim_::traits_type>
{
public:

    typedef Tsim_ simulator_type;
    typedef typename simulator_type::traits_type traits_type;
    typedef Domain<traits_type> base_type;
    typedef typename traits_type::world_type world_type;

    typedef typename world_type::particle_type particle_type;
    typedef typename world_type::particle_shape_type particle_shape_type;
    typedef typename world_type::molecule_info_type molecule_info_type;
    typedef typename world_type::species_id_type species_id_type;
    typedef typename world_type::position_type position_type;
    typedef typename world_type::particle_id_type particle_id_type;
    typedef typename world_type::length_type length_type;
    typedef typename world_type::size_type size_type;
    typedef typename world_type::structure_type structure_type;
    typedef typename world_type::particle_id_pair particle_id_pair;
    typedef typename world_type::particle_id_pair_and_distance
        particle_id_pair_and_distance;
    typedef typename world_type::particle_id_pair_and_distance_list
        particle_id_pair_and_distance_list;

    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::domain_id_type identifier_type;
    typedef typename traits_type::template shell_generator<
        ecell4::Sphere>::type spherical_shell_type;
    // typedef typename traits_type::template shell_generator<
    //     typename simulator_type::sphere_type>::type spherical_shell_type;
    typedef std::pair<const typename traits_type::shell_id_type, spherical_shell_type> spherical_shell_id_pair;
    typedef typename traits_type::reaction_record_type reaction_record_type;

    //XXX: typedef std::map<shell_id_type, spherical_shell_type> spherical_shell_map;
    typedef typename win32_workaround::get_mapper_mf<shell_id_type, spherical_shell_type>::type spherical_shell_map;
    typedef sized_iterator_range<typename spherical_shell_map::const_iterator> spherical_shell_id_pair_range;
    typedef MultiParticleContainer<traits_type> multi_particle_container_type;

    enum event_kind
    {
        NONE,
        ESCAPE,
        REACTION,
        NUM_MULTI_EVENT_KINDS
    };

private:
    struct last_reaction_setter: ReactionRecorder<reaction_record_type>
    {
        virtual ~last_reaction_setter() {}

        virtual void operator()(reaction_record_type const& rec)
        {
            outer_.last_reaction_.swap(const_cast<reaction_record_type&>(rec));
        }

        last_reaction_setter(Multi& outer): outer_(outer) {}

        Multi& outer_;
    };

    struct volume_clearer: VolumeClearer<particle_shape_type, particle_id_type>
    {
        virtual ~volume_clearer() {}

        virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore)
        {
            if (!outer_.within_shell(shape))
            {
                outer_.last_event_ = ESCAPE;
                return outer_.clear_volume(shape, ignore);
            }
            return true;
        }

        virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1)
        {
            if (!outer_.within_shell(shape))
            {
                outer_.last_event_ = ESCAPE;
                return outer_.clear_volume(shape, ignore0, ignore1);
            }
            return true;
        }

        volume_clearer(Multi& outer): outer_(outer) {}

        Multi& outer_;
    };

    friend struct volume_clearer;

public:
    virtual ~Multi() {}

    virtual char const* type_name() const
    {
        return "Multi";
    }

    virtual std::string as_string() const
    {
        return (boost::format(
            "%s(id=%s, event=%s, last_time=%.16g, dt=%.16g, particles=[%s])") %
            type_name() %
            boost::lexical_cast<std::string>(base_type::id_).c_str() %
            boost::lexical_cast<std::string>(base_type::event_.first).c_str() %
            base_type::last_time_ % base_type::dt_ %
            stringize_and_join(
                make_select_first_range(pc_.get_particles_range()),
                ", ")).str();
    }

    Multi(identifier_type const& id, simulator_type& main, Real dt_factor)
        : base_type(id), main_(main), pc_(*main.world()), dt_factor_(dt_factor),
          shells_(), last_event_(NONE)
    {
        BOOST_ASSERT(dt_factor > 0.);
        base_type::dt_ = dt_factor_ * determine_dt(*main_.world());
        // base_type::dt_ = dt_factor_ * BDSimulator<traits_type>::determine_dt(*main_.world());
    }

    static Real determine_dt(world_type const& world)
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(molecule_info_type s, world.get_molecule_info_range())
        {
            if (D_max < s.D)
                D_max = s.D;
            if (radius_min > s.radius)
                radius_min = s.radius;
        }
        return pow_2(radius_min * 2) / (D_max * 2);
    }

    event_kind const& last_event() const
    {
        return last_event_;
    }

    reaction_record_type const& last_reaction() const
    {
        return last_reaction_;
    }

    bool has_particle(particle_id_type const& pid) const
    {
        return pc_.has_particle(pid);
    }

    bool add_particle(particle_id_pair const& pp)
    {
        return pc_.update_particle(pp.first, pp.second);
    }

    bool add_shell(spherical_shell_id_pair const& sp)
    {
        spherical_shell_id_pair new_sp(sp);
        new_sp.second.did() = base_type::id();
        return shells_.insert(new_sp).second;
    }

    spherical_shell_id_pair_range get_shells() const
    {
        return spherical_shell_id_pair_range(shells_.begin(), shells_.end(), shells_.size());
    }

    virtual typename Domain<traits_type>::size_type num_shells() const
    {
        return shells_.size();
    }

    virtual typename Domain<traits_type>::size_type multiplicity() const
    {
        return pc_.num_particles();
    }

    virtual void accept(ImmutativeDomainVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeDomainVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    bool within_shell(particle_shape_type const& sphere) const
    {
        for (typename spherical_shell_map::const_iterator
                i(shells_.begin()), e(shells_.end()); i != e; ++i)
        {
            spherical_shell_id_pair const& sp(*i);
            position_type ppos(main_.world()->periodic_transpose(sphere.position(), (sp).second.position()));
            if (distance(ppos, (sp).second.shape().position()) < (sp).second.shape().radius() - sphere.radius())
            {
                return true;
            }
        }
        return false;
    }

    bool clear_volume(particle_shape_type const& shape, particle_id_type const& ignore) const
    {
        LOG_DEBUG(("clear_volume was called here."));
        main_.clear_volume(shape, base_type::id_);

        const particle_id_pair_and_distance_list overlapped(
            main_.world()->check_overlap(shape, ignore));
        if (overlapped.size() > 0)
        {
            return false;
        }
        return true;
        // return (main_.world()->no_overlap(shape, ignore));
    }

    bool clear_volume(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1) const
    {
        LOG_DEBUG(("clear_volume was called here."));
        main_.clear_volume(shape, base_type::id_);

        const particle_id_pair_and_distance_list overlapped(
            main_.world()->check_overlap(shape, ignore0, ignore1));
        if (overlapped.size() > 0)
        {
            return false;
        }
        return true;
        // return (main_.world()->no_overlap(shape, ignore0, ignore1));
    }

    typename multi_particle_container_type::particle_id_pair_range
    get_particles_range() const
    {
        return pc_.get_particles_range();
    }

    void step()
    {
        boost::scoped_ptr<
            typename multi_particle_container_type::transaction_type>
                tx(pc_.create_transaction());
        typedef typename multi_particle_container_type::transaction_type::particle_id_pair_generator particle_id_pair_generator;
        typedef typename multi_particle_container_type::transaction_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;
        last_reaction_setter rs(*this);
        volume_clearer vc(*this);
        BDPropagator<traits_type> ppg(
            *tx, *main_.network_rules(), main_.rng(),
            base_type::dt_,
            1 /* FIXME: dissociation_retry_moves */, &rs, &vc,
            make_select_first_range(pc_.get_particles_range()));

        last_event_ = NONE;

        while (ppg())
        {
            if (last_reaction_)
            {
                last_event_ = REACTION;
                break;
            }
        }
    }

protected:
    simulator_type& main_;
    multi_particle_container_type pc_;
    Real dt_factor_;
    spherical_shell_map shells_;
    event_kind last_event_;
    reaction_record_type last_reaction_;

    static Logger& log_;
};

template<typename Tsim_>
Logger& Multi<Tsim_>::log_(Logger::get_logger("ecell.Multi"));

#endif /* MULTI_HPP */
