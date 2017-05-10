#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>

#include "BDPropagator.hpp"
#include "World.hpp"
#include "ParticleSimulator.hpp"
#include "utils/pair.hpp"

#include "PotentialField.hpp"

#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/exceptions.hpp>

template<typename Tworld_>
struct BDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
};

template<typename Ttraits_>
class BDSimulator: public ParticleSimulator<Ttraits_>
{
public:

    typedef Ttraits_ traits_type;
    typedef ParticleSimulator<Ttraits_> base_type;

    typedef typename base_type::model_type model_type;

    typedef typename traits_type::world_type world_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::species_id_type species_id_type;
    // typedef typename world_type::species_type species_type;
    typedef typename world_type::molecule_info_type molecule_info_type;
    typedef typename world_type::particle_shape_type particle_shape_type;
    typedef typename world_type::particle_id_pair particle_id_pair;
    typedef typename world_type::particle_type particle_type;
    typedef typename world_type::particle_id_type particle_id_type;
    typedef typename world_type::traits_type::position_type position_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::rate_type rate_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename ReactionRecorderWrapper<reaction_record_type>::reaction_info_type reaction_info_type;

    typedef typename world_type::particle_container_type particle_container_type;
    typedef ecell4::PotentialField<particle_container_type> potential_field_type;
    typedef typename ecell4::utils::get_mapper_mf<species_id_type, boost::shared_ptr<potential_field_type> >::type potential_field_map_type;

public:

    Real const& dt_factor()
    {
        return dt_factor_;
    }

    virtual ~BDSimulator() {}

    BDSimulator(
        const boost::shared_ptr<world_type>& world,
        const boost::shared_ptr<model_type>& ecell4_model,
        Real bd_dt_factor = 1.0,
        int dissociation_retry_moves = 1)
        : base_type(world, ecell4_model),
          dt_factor_(bd_dt_factor),
          num_retries_(dissociation_retry_moves)
    {
        calculate_dt();
    }

    BDSimulator(
        const boost::shared_ptr<world_type>& world,
        Real bd_dt_factor = 1.0,
        int dissociation_retry_moves = 1)
        : base_type(world),
          dt_factor_(bd_dt_factor),
          num_retries_(dissociation_retry_moves)
    {
        calculate_dt();
    }

    virtual void initialize()
    {
        ;
    }

    virtual void calculate_dt()
    {
        set_dt(dt_factor_ * determine_dt());
        LOG_DEBUG(("dt=%f", base_type::dt_));
    }

    virtual void set_dt(const Real& dt)
    {
        base_type::dt_ = dt;
    }

    virtual void step()
    {
        _step(base_type::dt());
    }

    void add_potential(const ecell4::Species& sp, const Real& radius)
    {
        std::pair<typename potential_field_map_type::iterator, bool> retval
            = potentials_.insert(typename potential_field_map_type::value_type(sp, boost::shared_ptr<potential_field_type>(new ecell4::LeashPotentialField<typename potential_field_type::container_type>(radius))));
        if (!retval.second)
        {
            throw ecell4::AlreadyExists("never reach here.");
        }
    }

    void add_potential(const ecell4::Species& sp, const boost::shared_ptr<ecell4::Shape>& shape)
    {
        std::pair<typename potential_field_map_type::iterator, bool> retval
            = potentials_.insert(typename potential_field_map_type::value_type(sp, boost::shared_ptr<potential_field_type>(new ecell4::ShapedHardbodyPotentialField<typename potential_field_type::container_type>(shape))));
        if (!retval.second)
        {
            throw ecell4::AlreadyExists("never reach here.");
        }
    }

    void add_potential(const ecell4::Species& sp, const boost::shared_ptr<ecell4::Shape>& shape, const Real& threshold)
    {
        std::pair<typename potential_field_map_type::iterator, bool> retval
            = potentials_.insert(typename potential_field_map_type::value_type(sp, boost::shared_ptr<potential_field_type>(new ecell4::ShapedDiscretePotentialField<typename potential_field_type::container_type>(shape, threshold))));
        if (!retval.second)
        {
            throw ecell4::AlreadyExists("never reach here.");
        }
    }

    virtual bool step(const time_type& upto)
    {
        time_type const lt(upto - base_type::t());
        if (lt <= 0.0)
        {
            return false;
        }
        if (base_type::dt() < lt)
        {
            _step(base_type::dt());
        }
        else
        {
            _step(lt);
            base_type::set_t(upto);
        }
        return true;
    }

    Real determine_dt()
    {
        Real prob = 0.0;
        BOOST_FOREACH(reaction_rule_type const& rr,
                      (*base_type::network_rules_).zeroth_order_reaction_rules())
        {
            prob += rr.k();
        }
        prob *= (*base_type::world_).volume();

        if (prob == 0.0 && (*base_type::world_).num_particles() == 0)
        {
            return std::numeric_limits<Real>::infinity();
        }

        Real D_max(0.0), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(molecule_info_type info,
                      (*base_type::world_).get_molecule_info_range())
        {
            if (D_max < info.D)
            {
                D_max = info.D;
            }
            if (radius_min > info.radius)
            {
                radius_min = info.radius;
            }
        }
        const Real dt(gsl_pow_2(radius_min * 2) / (D_max * 2));
        return (prob == 0.0 ? dt : std::min(dt, 1.0 / prob));
    }

    virtual bool check_reaction() const
    {
        return last_reactions().size() > 0;
    }

    std::vector<std::pair<ecell4::ReactionRule, reaction_info_type> > last_reactions() const
    {
        return (*dynamic_cast<ReactionRecorderWrapper<reaction_record_type>*>(
            base_type::rrec_.get())).last_reactions();
    }

protected:

    void _step(time_type dt)
    {
        {
            BDPropagator<traits_type> propagator(
                *base_type::world_,
                *base_type::network_rules_,
                base_type::rng(),
                dt, num_retries_,
                base_type::rrec_.get(), 0,
                make_select_first_range(base_type::world_->
                                        get_particles_range()),
                potentials_);
            while (propagator());
        }

        try
        {
            attempt_zeroth_order_reaction(dt);
        }
        catch (no_space const&)
        {
            LOG_DEBUG(("birth reaction rejected."));
            // ++rejected_moves_;
        }

        LOG_DEBUG(("%d: t=%lg, dt=%lg", base_type::num_steps_,
                   base_type::t(), dt));

        ++base_type::num_steps_;
        base_type::set_t(base_type::t() + dt);
    }

    bool attempt_zeroth_order_reaction(time_type const dt)
    {
        typename network_rules_type::reaction_rules const
            rules((*base_type::network_rules_).zeroth_order_reaction_rules());
        if (::size(rules) == 0)
            return false;

        const Real rnd(
            base_type::rng().random() / (dt * (*base_type::world_).volume()));
        Real prob = 0.0;

        BOOST_FOREACH(reaction_rule_type const& rr, rules)
        {
            prob += rr.k();
            if (prob > rnd)
            {
                typename reaction_rule_type::species_id_range products(
                        rr.get_products());
                BOOST_ASSERT(::size(products) == 1);
                const molecule_info_type minfo(
                    (*base_type::world_).get_molecule_info(products[0]));

                //XXX: A cuboidal region is expected here.
                const position_type& edge_lengths((*base_type::world_).edge_lengths());
                const position_type new_pos(
                    base_type::rng().uniform(0, edge_lengths[0]),
                    base_type::rng().uniform(0, edge_lengths[1]),
                    base_type::rng().uniform(0, edge_lengths[2]));

                const particle_shape_type new_particle(new_pos, minfo.radius);
                if (!(*base_type::world_).no_overlap(new_particle))
                {
                    LOG_INFO(("no space for product particle."));
                    throw no_space();
                }

                particle_id_pair pp(
                    (*base_type::world_).new_particle(products[0], new_pos).first);

                if (base_type::rrec_)
                {
                    (*base_type::rrec_)(
                        reaction_record_type(rr.id(), array_gen(pp)));
                }
                return true;
            }
        }
        return false;
    }

private:

    Real const dt_factor_;
    int const num_retries_;
    Real R_;
    static Logger& log_;

    potential_field_map_type potentials_;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));

#endif /* BD_SIMULATOR_HPP */
