#ifndef MULTI_HPP
#define MULTI_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/algorithm/string/join.hpp>

#include "exceptions.hpp"
#include "Domain.hpp"
#include "ParticleContainer.hpp"
#include "ParticleContainerBase.hpp"
#include "Sphere.hpp"
#include "BDSimulator.hpp"
#include "BDPropagator.hpp"
#include "Logger.hpp"
#include "PairGreensFunction.hpp"
#include "utils/array_helper.hpp"
#include "utils/range.hpp"

template<typename Ttraits_>
class MultiParticleContainer
    : public Ttraits_::world_type::particle_container_type
{
public:
    typedef ParticleContainerUtils<typename Ttraits_::world_type::traits_type> utils;
    typedef typename Ttraits_::world_type world_type;
    typedef typename world_type::traits_type traits_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename particle_type::shape_type particle_shape_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef Transaction<traits_type> transaction_type;
    typedef abstract_limited_generator<particle_id_pair> particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;
    typedef std::map<particle_id_type, particle_type> particle_map;
    typedef sized_iterator_range<typename particle_map::const_iterator> particle_id_pair_range;

    virtual ~MultiParticleContainer() {}

    virtual size_type num_particles() const
    {
        return particles_.size();
    }

    virtual length_type world_size() const
    {
        return world_.world_size();
    }

    virtual species_type const& get_species(species_id_type const& id) const
    {
        return world_.get_species(id);
    }

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return world_.get_structure(id);
    }

    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        particle_id_pair const retval(world_.new_particle(sid, pos));
        particles_.insert(retval);
        return retval;
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        bool const retval(world_.update_particle(pi_pair));
        particles_[pi_pair.first] = pi_pair.second;
        return retval;
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        world_.remove_particle(id);
        return particles_.erase(id);
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

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return check_overlap(s, array_gen<particle_id_type>());
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return check_overlap(s, array_gen(ignore));
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return check_overlap(s, array_gen(ignore1, ignore2));
    }

    template<typename Tsph_, typename Tset_>
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s, Tset_ const& ignore) const
    {
        typename utils::template overlap_checker<Tset_> checker(ignore);
        for (typename particle_map::const_iterator i(particles_.begin()),
                                                   e(particles_.end());
             i != e; ++i)
        {
            length_type const dist(world_.distance(shape((*i).second), s.position()));
            if (dist < s.radius())
            {
                checker(i, dist);
            }
        }
        return checker.result();
    }

    virtual particle_id_pair_generator* get_particles() const
    {
        return make_range_generator<particle_id_pair>(particles_);
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

    virtual length_type apply_boundary(length_type const& v) const
    {
        return world_.apply_boundary(v);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return world_.cyclic_transpose(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return world_.cyclic_transpose(p0, p1);
    }

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range(particles_.begin(), particles_.end(),
                                      particles_.size());
    }

    MultiParticleContainer(world_type& world): world_(world) {}

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
    typedef typename traits_type::world_type::particle_type particle_type;
    typedef typename particle_type::shape_type particle_shape_type;
    typedef typename traits_type::world_type::species_type species_type;
    typedef typename traits_type::world_type::species_id_type species_id_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_type particle_id_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::size_type size_type;
    typedef typename traits_type::world_type::structure_type structure_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::domain_id_type identifier_type;
    typedef typename traits_type::template shell_generator<
        typename simulator_type::sphere_type>::type spherical_shell_type;
    typedef std::pair<const typename traits_type::shell_id_type, spherical_shell_type> spherical_shell_id_pair;
    typedef typename traits_type::reaction_record_type reaction_record_type;

    typedef std::map<shell_id_type, spherical_shell_type> spherical_shell_map;
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

public:
    virtual ~Multi() {}

    virtual char const* type_name() const
    {
        return "Multi";
    }

    virtual std::string as_string() const
    {
        return (boost::format(
            "%s(id=%s, event=%s, last_time=%g, dt=%g, particles=[%s])") %
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
          shells_(), last_event_(NONE) {}

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
        return pc_.update_particle(pp);
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
            position_type ppos(main_.world()->cyclic_transpose(sphere.position(), (sp).second.position()));
            if (distance(ppos, (sp).second.shape().position()) < (sp).second.shape().radius())
            {
                return true;
            }
        }
        return false;
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
        BDPropagator<traits_type> ppg(
            *tx, *main_.network_rules(), main_.rng(),
            dt_factor_ * BDSimulator<traits_type>::determine_dt(*main_.world()),
            1 /* FIXME: dissociation_retry_moves */, &rs,
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

        boost::scoped_ptr<particle_id_pair_generator>
            added_particles(tx->get_added_particles()),
            modified_particles(tx->get_modified_particles());
        LOG_DEBUG(("added_particles=%zu, modified_particles=%zu",
            ::count(*added_particles), ::count(*modified_particles)));
        chained_generator<particle_id_pair_generator,
                          particle_id_pair_generator>
            gen(*added_particles, *modified_particles);
        while (::valid(gen))
        {
            particle_id_pair pp(gen());
            boost::scoped_ptr<particle_id_pair_and_distance_list>
                overlapped(main_.world()->check_overlap(
                    pp.second.shape(), pp.first));
            if (overlapped && ::size(*overlapped))
            {
                log_.info("collision occurred between particles of a multi and the outside: %s.  moves will be rolled back.",
                    boost::lexical_cast<std::string>(pp.first).c_str());
                tx->rollback();
                return;
            }

            if (!within_shell(pp.second.shape()))
            {
                last_event_ = ESCAPE;
                main_.clear_volume(pp.second.shape(), base_type::id_);
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
