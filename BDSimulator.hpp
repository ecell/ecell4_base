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
    typedef typename traits_type::random_number_engine random_number_engine;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::world_type::species_id_type species_id_type;
    typedef typename traits_type::world_type::species_type species_type;
    typedef typename traits_type::rate_type rate_type;
    typedef typename traits_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::network_rules_type network_rules_type;

public:
    Real const& dt_factor()
    {
        return dt_factor_;
    }

    virtual ~BDSimulator() {}

    BDSimulator(world_type& world, rng_type& rng,
                network_rules_type const& network_rules)
        : base_type(world, rng, network_rules),
          num_steps_(0), num_reactions_(0), dt_factor(1e-5) {}

    virtual void initialize()
    {
        determine_dt();
        LOG_DEBUG(("dt=%f", dt_));
    }

    virtual void step()
    {
        {
            boost::scoped_ptr<typename world_type::transaction_type> tx(
                    world_.create_transaction());
            BDPropagator<traits_type> propagator(world_, *tx, *this);
            while (propagator());
            num_reactions_ += propagator.get_reactions().size();
            LOG_DEBUG(("%d: t=%lg, dt=%lg, reactions=%d", num_steps_, t_, dt_, num_reactions_));
        }
        ++num_steps_;
        t_ += dt_;
    }

protected:
    void determine_dt()
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, world_.get_species())
        {
            if (D_max < s.D())
                D_max = s.D();
            if (radius_min > s.radius())
                radius_min = s.radius();
        }
        dt_ = dt_factor_ * gsl_pow_2(radius_min * 2) / D_max * 2;
    }

private:
    int num_steps_;
    int num_reactions_;
    Real dt_factor_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));


#endif /* BD_SIMULATOR_HPP */
