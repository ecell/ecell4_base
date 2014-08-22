#ifndef __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP
#define __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "GillespieWorld.hpp"


namespace ecell4
{

namespace gillespie
{

class ReactionRuleEvent
{
public:

    ReactionRuleEvent()
        : rr_()
    {
        ;
    }

    ReactionRuleEvent(const ReactionRule& rr)
        : rr_(rr)
    {
        ;
    }

    virtual ~ReactionRuleEvent()
    {
        ;
    }

    inline const Integer get_coef(const Species& pttrn, const Species& sp) const
    {
        return SpeciesExpressionMatcher(pttrn).count(sp);
    }

    virtual void initialize(GillespieWorld* w) = 0;
    virtual void inc(const Species& sp, const Integer val = +1) = 0;
    virtual const Real propensity(const Real& volume) const = 0;

    inline void dec(const Species& sp)
    {
        inc(sp, -1);
    }

    ReactionRule draw(GillespieWorld* w)
    {
        const std::pair<ReactionRule::reactant_container_type, Integer>
            retval(__draw(w));
        const std::vector<std::vector<Species> >
            possibles(rrgenerate(rr_, retval.first));

        assert(retval.second > 0);
        assert(retval.second >= possibles.size());

        if (possibles.size() == 0)
        {
            return ReactionRule();
        }
        else if (retval.second == 1)
        {
            // assert(possibles.size() == 1);
            return ReactionRule(retval.first, possibles[0], rr_.k());
        }
        else
        {
            const Integer rnd2(w->rng()->uniform_int(0, retval.second - 1));
            if (rnd2 >= possibles.size())
            {
                return ReactionRule();
            }
            return ReactionRule(retval.first, possibles[rnd2], rr_.k());
        }
    }

protected:

    virtual std::pair<ReactionRule::reactant_container_type, Integer>
        __draw(GillespieWorld* w) = 0;

protected:

    ReactionRule rr_;
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

    FirstOrderReactionRuleEvent(const ReactionRule& rr)
        : base_type(rr), num_tot1_(0)
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

    void initialize(GillespieWorld* w)
    {
        const std::vector<Species>& species(w->list_species());
        const ReactionRule::reactant_container_type& reactants(rr_.reactants());

        num_tot1_ = 0;
        for (std::vector<Species>::const_iterator i(species.begin());
            i != species.end(); ++i)
        {
            const Integer coef(get_coef(reactants[0], *i));
            if (coef > 0)
            {
                num_tot1_ += coef * w->num_molecules_exact(*i);
            }
        }
    }

    std::pair<ReactionRule::reactant_container_type, Integer>
    __draw(GillespieWorld* w)
    {
        const std::vector<Species>& species(w->list_species());
        const ReactionRule::reactant_container_type& reactants(rr_.reactants());

        const Real rnd1(w->rng()->uniform(0.0, num_tot1_));

        Integer num_tot(0);
        for (std::vector<Species>::const_iterator i(species.begin());
            i != species.end(); ++i)
        {
            const Integer coef(get_coef(reactants[0], *i));
            if (coef > 0)
            {
                num_tot += coef * w->num_molecules_exact(*i);
                if (num_tot >= rnd1)
                {
                    return std::make_pair(
                        ReactionRule::reactant_container_type(1, *i), coef);
                }
            }
        }

        ; // never get here
    }

    const Real propensity(const Real& volume) const
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

    SecondOrderReactionRuleEvent(const ReactionRule& rr)
        : base_type(rr), num_tot1_(0), num_tot2_(0), num_tot12_(0)
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

    void initialize(GillespieWorld* w)
    {
        const std::vector<Species>& species(w->list_species());
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
                const Integer num(w->num_molecules_exact(*i));
                const Integer tmp(coef1 * num);
                num_tot1_ += tmp;
                num_tot2_ += coef2 * num;
                num_tot12_ += coef2 * tmp;
            }
        }
    }

    std::pair<ReactionRule::reactant_container_type, Integer>
    __draw(GillespieWorld* w)
    {
        const std::vector<Species>& species(w->list_species());
        const ReactionRule::reactant_container_type& reactants(rr_.reactants());

        const Real rnd1(w->rng()->uniform(0.0, num_tot1_));

        Integer num_tot(0), coef1(0);
        std::vector<Species>::const_iterator itr1(species.begin());
        for (; itr1 != species.end(); ++itr1)
        {
            const Integer coef(get_coef(reactants[0], *itr1));
            if (coef > 0)
            {
                num_tot += coef * w->num_molecules_exact(*itr1);
                if (num_tot >= rnd1)
                {
                    coef1 = coef;
                    break;
                }
            }
        }

        const Real rnd2(
            w->rng()->uniform(0.0, num_tot2_ - get_coef(reactants[0], *itr1)));

        num_tot = 0;
        for (std::vector<Species>::const_iterator i(species.begin());
            i != species.end(); ++i)
        {
            const Integer coef(get_coef(reactants[1], *i));
            if (coef > 0)
            {
                const Integer num(w->num_molecules_exact(*i));
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

        ; // never get here
    }

    const Real propensity(const Real& volume) const
    {
        return (num_tot1_ * num_tot2_ - num_tot12_) * rr_.k() / volume;
    }

protected:

    Integer num_tot1_, num_tot2_, num_tot12_;
};

class GillespieSimulator
    : public Simulator<Model, GillespieWorld>
{
public:

    typedef Simulator<Model, GillespieWorld> base_type;

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

    Real t(void) const;
    Real dt(void) const;

    void step(void) ;
    bool step(const Real & upto);

    // Optional members

    void set_t(const Real &t);
    std::vector<ReactionRule> last_reactions() const;

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
    // Integer num_molecules(
    //     const Model::reaction_rule_container_type::size_type& u);
    // ReactionRule draw_exact_reaction(
    //     const Model::reaction_rule_container_type::size_type& u);
    // std::pair<ReactionRule::reactant_container_type, Integer>
    //     draw_exact_reactants(const Model::reaction_rule_container_type::size_type& u);
    // std::pair<ReactionRule::reactant_container_type, Integer>
    //     draw_exact_reactants(const Species& sp1, const stoichiometry_container_type& retval);
    // std::pair<ReactionRule::reactant_container_type, Integer>
    //     draw_exact_reactants(const Species& sp1, const Species& sp2, const stoichiometry_container_type& retval);

    // Real calculate_propensity(
    //     const Model::reaction_rule_container_type::size_type& u);
    // stoichiometry_container_type get_stoichiometry(const Species& sp);
    // stoichiometry_container_type get_stoichiometry(
    //     const Species& sp1, const Species& sp2);

    // void calculate_stoichiometries();
    // void append_stoichiometries(const Species& sp);
    void add_molecules(const Species& sp);
    void remove_molecules(const Species& sp);

protected:

    Real dt_;
    ReactionRule next_reaction_;
    std::vector<ReactionRule> last_reactions_;

    boost::ptr_vector<ReactionRuleEvent> events_;
};

}

} // ecell4

#endif /* __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP */
