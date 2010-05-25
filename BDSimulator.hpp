#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "NetworkRules.hpp"
#include "BDPropagator.hpp"
#include "World.hpp"
#include "ParticleSimulator.hpp"
#include "utils/pair.hpp"

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
    typedef typename traits_type::world_type world_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::species_id_type species_id_type;
    typedef typename world_type::species_type species_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::rate_type rate_type;

public:
    Real const& dt_factor()
    {
        return dt_factor_;
    }

    virtual ~BDSimulator() {}

    BDSimulator(world_type& world, 
                network_rules_type const& network_rules,
                rng_type& rng,
                Real dt_factor = .5,
                int dissociation_retry_moves = 1)
        : base_type(world, network_rules, rng),
          dt_factor_(dt_factor), num_retrys_(dissociation_retry_moves)
    {
        determine_dt();
        LOG_DEBUG(("dt=%f", dt_));
    }

    virtual void step()
    {
        step(base_type::dt_);
    }

    virtual bool step(time_type const& upto)
    {
        time_type const lt(upto - base_type::t_);
        if (lt <= 0.)
            return false;
        if (base_type::dt_ < lt)
        {
            _step(base_type::dt_);
        }
        else
        {
            _step(lt);
            base_type::t_ = upto;
        }
        return true;
    }

protected:
    void determine_dt()
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, base_type::world_.get_species())
        {
            if (D_max < s.D())
                D_max = s.D();
            if (radius_min > s.radius())
                radius_min = s.radius();
        }
        base_type::dt_ = dt_factor_ * gsl_pow_2(radius_min * 2) / (D_max * 2);
    }

    void _step(time_type const& dt)
    {
        {
            BDPropagator<traits_type> propagator(
                base_type::world_,
                base_type::network_rules_,
                base_type::rng_,
                dt, num_retrys_,
                make_select_first_range(base_type::world_.get_particles_range()));
            while (propagator());
            base_type::num_reactions_ += propagator.get_reactions().size();
            LOG_DEBUG(("%d: t=%lg, dt=%lg, reactions=%d", base_type::num_steps_, t_, dt, base_type::num_reactions_));
        }
        ++base_type::num_steps_;
        base_type::t_ += dt;
    }

private:
    Real const dt_factor_;
    int const num_retrys_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));


#endif /* BD_SIMULATOR_HPP */
