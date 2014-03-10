#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>
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
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type 
    reaction_recorder_type;

public:
    Real const& dt_factor()
    {
        return dt_factor_;
    }

    virtual ~BDSimulator() {}

    BDSimulator(boost::shared_ptr<world_type> world, 
                boost::shared_ptr<network_rules_type const> network_rules,
                rng_type& rng, Real dt_factor = 1.,
                int dissociation_retry_moves = 1)
        : base_type(world, network_rules, rng),
          dt_factor_(dt_factor), num_retries_(dissociation_retry_moves)
    {
        calculate_dt();
    }

    virtual void calculate_dt()
    {
        base_type::dt_ = dt_factor_ * determine_dt(*base_type::world_);
        LOG_DEBUG(("dt=%f", base_type::dt_));
    }

    virtual void step()
    {
        _step(base_type::dt_);
    }

    virtual bool step(time_type upto)
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

    static Real determine_dt(world_type const& world)
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, world.get_species())
        {
            if (D_max < s.D())
                D_max = s.D();
            if (radius_min > s.radius())
                radius_min = s.radius();
        }
        return gsl_pow_2(radius_min * 2) / (D_max * 2);
    }

protected:
    void _step(time_type dt)
    {
        {
            BDPropagator<traits_type> propagator(
                *base_type::world_,
                *base_type::network_rules_,
                base_type::rng_,
                dt, num_retries_,
                base_type::rrec_.get(), 0,
                make_select_first_range(base_type::world_->
                                        get_particles_range()));
            while (propagator());
            LOG_DEBUG(("%d: t=%lg, dt=%lg", base_type::num_steps_, 
                       base_type::t_, dt));
        }
        ++base_type::num_steps_;
        base_type::t_ += dt;
    }

private:
    Real const dt_factor_;
    int const num_retries_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));


#endif /* BD_SIMULATOR_HPP */
