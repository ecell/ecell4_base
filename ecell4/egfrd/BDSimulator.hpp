#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>

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

    typedef typename base_type::model_type model_type;

    typedef typename traits_type::world_type world_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::species_id_type species_id_type;
    // typedef typename world_type::species_type species_type;
    typedef typename world_type::molecule_info_type molecule_info_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::rate_type rate_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;

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

    virtual void initialize()
    {
        ;
    }

    virtual void calculate_dt()
    {
        set_dt(dt_factor_ * determine_dt(*base_type::world_));
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

    Real determine_dt(world_type const& world)
    {
        Real D_max(0.0), radius_min(std::numeric_limits<Real>::max());

        // BOOST_FOREACH(species_type s, world.get_species())
        // {
        //     if (D_max < s.D())
        //         D_max = s.D();
        //     if (radius_min > s.radius())
        //         radius_min = s.radius();
        // }
        BOOST_FOREACH(molecule_info_type info, (*base_type::world_).get_molecule_info_range())
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
        return gsl_pow_2(radius_min * 2) / (D_max * 2);
    }

    virtual bool check_reaction() const
    {
        return last_reactions().size() > 0;
    }

    std::vector<ecell4::ReactionRule> last_reactions() const
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
                                        get_particles_range()));
            while (propagator());
            LOG_DEBUG(("%d: t=%lg, dt=%lg", base_type::num_steps_,
                       base_type::t(), dt));
        }
        ++base_type::num_steps_;
        base_type::set_t(base_type::t() + dt);
    }

private:

    Real const dt_factor_;
    int const num_retries_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));

#endif /* BD_SIMULATOR_HPP */
