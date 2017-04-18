#ifndef ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP
#define ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/SimulatorBase.hpp>

#include "GillespieWorld.hpp"


namespace ecell4
{

namespace gillespie
{

class ReactionInfo
{
public:

    typedef std::vector<Species> container_type;

public:

    ReactionInfo(
        const Real t,
        const container_type& reactants,
        const container_type& products)
        : t_(t), reactants_(reactants), products_(products)
    {}

    ReactionInfo(
        const ReactionInfo& another)
        : t_(another.t()), reactants_(another.reactants()), products_(another.products())
    {}

    Real t() const
    {
        return t_;
    }

    const container_type& reactants() const
    {
        return reactants_;
    }

    void add_reactant(const Species& sp)
    {
        reactants_.push_back(sp);
    }

    const container_type& products() const
    {
        return products_;
    }

    void add_product(const Species& sp)
    {
        products_.push_back(sp);
    }

protected:

    Real t_;
    container_type reactants_, products_;
};

class GillespieSimulator
    : public SimulatorBase<Model, GillespieWorld>
{
public:

    typedef SimulatorBase<Model, GillespieWorld> base_type;
    typedef ReactionInfo reaction_info_type;

protected:

    class ReactionRuleEvent
    {
    public:

        ReactionRuleEvent()
            : sim_(), rr_()
        {
            ;
        }

        ReactionRuleEvent(GillespieSimulator* sim, const ReactionRule& rr)
            : sim_(sim), rr_(rr)
        {
            ;
        }

        virtual ~ReactionRuleEvent()
        {
            ;
        }

        const ReactionRule& reaction_rule() const
        {
            return rr_;
        }

        inline const Integer get_coef(const Species& pttrn, const Species& sp) const
        {
            return sim_->model()->apply(pttrn, sp);
        }

        inline const std::vector<ReactionRule> generate(
            const ReactionRule::reactant_container_type& reactants) const
        {
            return sim_->model()->apply(rr_, reactants);
        }

        virtual void initialize() = 0;
        virtual void inc(const Species& sp, const Integer val = +1) = 0;
        virtual const Real propensity() const = 0;

        inline void dec(const Species& sp)
        {
            inc(sp, -1);
        }

        ReactionRule draw()
        {
            const std::pair<ReactionRule::reactant_container_type, Integer>
                retval(__draw());
            if (retval.second == 0)
            {
                return ReactionRule();
            }

            const std::vector<ReactionRule> reactions(generate(retval.first));

            assert(retval.second > 0);
            assert(retval.second >= reactions.size());

            if (reactions.size() == 0)
            {
                return ReactionRule();
            }
            else if (retval.second == 1)
            {
                // assert(possibles.size() == 1);
                return reactions[0];
            }
            else
            {
                const ReactionRule::reactant_container_type::size_type rnd2(
                    static_cast<ReactionRule::reactant_container_type::size_type>(
                        rng()->uniform_int(0, retval.second - 1)));
                if (rnd2 >= reactions.size())
                {
                    return ReactionRule();
                }
                return reactions[rnd2];
            }
        }

    protected:

        inline const boost::shared_ptr<RandomNumberGenerator>& rng() const
        {
            return sim_->world()->rng();
        }

        inline const GillespieWorld& world() const
        {
            return (*sim_->world());
        }

        virtual std::pair<ReactionRule::reactant_container_type, Integer>
            __draw() = 0;

    protected:

        GillespieSimulator* sim_;
        ReactionRule rr_;
    };

    class ZerothOrderReactionRuleEvent
        : public ReactionRuleEvent
    {
    public:

        typedef ReactionRuleEvent base_type;

        ZerothOrderReactionRuleEvent()
            : base_type()
        {
            ;
        }

        ZerothOrderReactionRuleEvent(GillespieSimulator* sim, const ReactionRule& rr)
            : base_type(sim, rr)
        {
            ;
        }

        void inc(const Species& sp, const Integer val = +1)
        {
            ; // do nothing
        }

        void initialize()
        {
            ; // do nothing
        }

        std::pair<ReactionRule::reactant_container_type, Integer> __draw()
        {
            return std::make_pair(ReactionRule::reactant_container_type(), 1);
        }

        const Real propensity() const
        {
            return rr_.k() * sim_->world()->volume();
        }
    };


    class FirstOrderReactionRuleEvent
        : public ReactionRuleEvent
    {
    public:

        typedef ReactionRuleEvent base_type;

        FirstOrderReactionRuleEvent()
            : base_type(), num_tot1_(0)
        {
            ;
        }

        FirstOrderReactionRuleEvent(GillespieSimulator* sim, const ReactionRule& rr)
            : base_type(sim, rr), num_tot1_(0)
        {
            ;
        }

        void inc(const Species& sp, const Integer val = +1)
        {
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());
            const Integer coef(get_coef(reactants[0], sp));
            if (coef > 0)
            {
                num_tot1_ += coef * val;
            }
        }

        void initialize()
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            num_tot1_ = 0;
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[0], *i));
                if (coef > 0)
                {
                    num_tot1_ += coef * world().num_molecules_exact(*i);
                }
            }
        }

        std::pair<ReactionRule::reactant_container_type, Integer> __draw()
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            const Real rnd1(rng()->uniform(0.0, num_tot1_));

            Integer num_tot(0);
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[0], *i));
                if (coef > 0)
                {
                    num_tot += coef * world().num_molecules_exact(*i);
                    if (num_tot >= rnd1)
                    {
                        return std::make_pair(
                            ReactionRule::reactant_container_type(1, *i), coef);
                    }
                }
            }

            return std::make_pair(ReactionRule::reactant_container_type(), 0);
        }

        const Real propensity() const
        {
            return num_tot1_ * rr_.k();
        }

    protected:

        Integer num_tot1_;
    };

    class SecondOrderReactionRuleEvent:
        public ReactionRuleEvent
    {
    public:

        typedef ReactionRuleEvent base_type;

        SecondOrderReactionRuleEvent()
            : base_type(), num_tot1_(0), num_tot2_(0), num_tot12_(0)
        {
            ;
        }

        SecondOrderReactionRuleEvent(GillespieSimulator* sim, const ReactionRule& rr)
            : base_type(sim, rr), num_tot1_(0), num_tot2_(0), num_tot12_(0)
        {
            ;
        }

        void inc(const Species& sp, const Integer val = +1)
        {
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());
            const Integer coef1(get_coef(reactants[0], sp));
            const Integer coef2(get_coef(reactants[1], sp));
            if (coef1 > 0 || coef2 > 0)
            {
                const Integer tmp(coef1 * val);
                num_tot1_ += tmp;
                num_tot2_ += coef2 * val;
                num_tot12_ += coef2 * tmp;
            }
        }

        void initialize()
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            num_tot1_ = 0;
            num_tot2_ = 0;
            num_tot12_ = 0;
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef1(get_coef(reactants[0], *i));
                const Integer coef2(get_coef(reactants[1], *i));
                if (coef1 > 0 || coef2 > 0)
                {
                    const Integer num(world().num_molecules_exact(*i));
                    const Integer tmp(coef1 * num);
                    num_tot1_ += tmp;
                    num_tot2_ += coef2 * num;
                    num_tot12_ += coef2 * tmp;
                }
            }
        }

        std::pair<ReactionRule::reactant_container_type, Integer> __draw()
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            const Real rnd1(rng()->uniform(0.0, num_tot1_));

            Integer num_tot(0), coef1(0);
            std::vector<Species>::const_iterator itr1(species.begin());
            for (; itr1 != species.end(); ++itr1)
            {
                const Integer coef(get_coef(reactants[0], *itr1));
                if (coef > 0)
                {
                    num_tot += coef * world().num_molecules_exact(*itr1);
                    if (num_tot >= rnd1)
                    {
                        coef1 = coef;
                        break;
                    }
                }
            }

            const Real rnd2(
                rng()->uniform(0.0, num_tot2_ - get_coef(reactants[0], *itr1)));

            num_tot = 0;
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[1], *i));
                if (coef > 0)
                {
                    const Integer num(world().num_molecules_exact(*i));
                    num_tot += coef * (i == itr1 ? num - 1 : num);
                    if (num_tot >= rnd2)
                    {
                        ReactionRule::reactant_container_type exact_reactants(2);
                        exact_reactants[0] = *itr1;
                        exact_reactants[1] = *i;
                        return std::make_pair(exact_reactants, coef1 * coef);
                    }
                }
            }

            return std::make_pair(ReactionRule::reactant_container_type(), 0);
        }

        const Real propensity() const
        {
            return (num_tot1_ * num_tot2_ - num_tot12_) * rr_.k() / world().volume();
        }

    protected:

        Integer num_tot1_, num_tot2_, num_tot12_;
    };

public:

    GillespieSimulator(
        boost::shared_ptr<Model> model,
        boost::shared_ptr<GillespieWorld> world)
        : base_type(model, world)
    {
        initialize();
    }

    GillespieSimulator(boost::shared_ptr<GillespieWorld> world)
        : base_type(world)
    {
        initialize();
    }

    // SimulatorTraits
    Real dt(void) const;

    void step(void) ;
    bool step(const Real & upto);

    // Optional members

    virtual bool check_reaction() const
    {
        return last_reactions_.size() > 0;
    }

    std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions() const
    {
        return last_reactions_;
    }

    /**
     * recalculate reaction propensities and draw the next time.
     */
    void initialize();

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    bool __draw_next_reaction(void);
    void draw_next_reaction(void);
    void increment_molecules(const Species& sp);
    void decrement_molecules(const Species& sp);

protected:

    Real dt_;
    ReactionRule next_reaction_rule_, next_reaction_;
    std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions_;

    boost::ptr_vector<ReactionRuleEvent> events_;
};

}

} // ecell4

#endif /* ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP */
