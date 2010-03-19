#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "NetworkRules.hpp"
#include "BDPropagator.hpp"
#include "Model.hpp"
#include "World.hpp"
#include "ParticleSimulator.hpp"

template<typename Tworld_>
struct BDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
};

template<typename Ttraits_>
class BDSimulator
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type world_type;
    typedef typename traits_type::model_type model_type;
    typedef typename traits_type::random_number_engine random_number_engine;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::world_type::species_id_type species_id_type;
    typedef typename traits_type::world_type::species_type species_type;
    typedef typename traits_type::rate_type rate_type;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef ReactionRuleInfo<
        typename reaction_rule_type::identifier_type,
        species_id_type, rate_type> reaction_rule_info_type;
    typedef NetworkRulesWrapper<network_rules_type, reaction_rule_info_type> network_rules_wrapper_type;

public:
    static const Real default_dt_factor = 1e-5;

public:
    model_type const& get_model() const
    {
        return model_;
    }

    world_type& get_world()
    {
        BOOST_ASSERT(world_);
        return *world_;
    }

    world_type const& get_world() const
    {
        BOOST_ASSERT(world_);
        return *world_;
    }

    network_rules_type const& get_network_rules() const
    {
        return model_.network_rules();
    }

    random_number_engine& get_rng()
    {
        return rng_;
    }

    Real get_t() const
    {
        return t_;
    }

    void set_t(Real val)
    {
        t_ = val;
    }

    Real get_dt() const
    {
        return dt_;
    }

    void set_dt(Real val)
    {
        dt_ = val;
    }

    Real get_dt_factor()
    {
        return dt_factor_;
    }

    void set_dt_factor(Real val)
    {
        dt_factor_ = val;
    }

    ~BDSimulator()
    {
        delete world_;
    }

    BDSimulator(Model const& model): model_(model), num_steps_(0), num_reactions_(0), world_(0)
    {
        set_dt_factor(1e-5);
        set_dt(0.);
        set_t(0.);
        setup();
    }

    void initialize()
    {
        determine_dt();
        LOG_DEBUG(("dt=%f", get_dt()));
    }

    void step()
    {
        {
            boost::scoped_ptr<typename world_type::transaction_type> tx(
                    get_world().create_transaction());
            BDPropagator<traits_type> propagator(*world_, *this, *tx);
            while (propagator());
            num_reactions_ += propagator.get_reactions().size();
            LOG_DEBUG(("%d: t=%lg, dt=%lg, reactions=%d", num_steps_, get_t(), get_dt(), num_reactions_));
        }
        ++num_steps_;
        set_t(get_t() + get_dt());
    }

protected:
    void determine_dt()
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, get_world().get_species())
        {
            if (D_max < s.D())
                D_max = s.D();
            if (radius_min > s.radius())
                radius_min = s.radius();
        }
        set_dt(get_dt_factor() * gsl_pow_2(radius_min * 2) / D_max * 2);
    }

    void setup()
    {
        typename world_type::length_type world_size(1.0);
        typename world_type::size_type matrix_size(10);

        try
        {
            world_size = boost::lexical_cast<typename world_type::length_type>(model_["world_size"]);
        }
        catch (...) {}

        try
        {
            matrix_size = boost::lexical_cast<typename world_type::size_type>(model_["matrix_size"]);
        }
        catch (...) {}

        world_ = new world_type(world_size, matrix_size);

        BOOST_FOREACH(typename model_type::species_type_type* st,
                      model_.get_species_types())
        {
            world_->add_species(
                species_type(st->id(),
                    boost::lexical_cast<Real>((*st)["D"]), 
                    boost::lexical_cast<Real>((*st)["radius"])));
        }
    }

private:
    Model const& model_;
    Real t_;
    Real dt_;
    Real dt_factor_;
    random_number_engine rng_;
    int num_steps_;
    int num_reactions_;
    world_type* world_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));


#endif /* BD_SIMULATOR_HPP */
