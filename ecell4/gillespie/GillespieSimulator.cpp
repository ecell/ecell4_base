#include "GillespieSimulator.hpp"
#include <numeric>
#include <vector>
#include <gsl/gsl_sf_log.h>

#include <cstring>
#include <sstream>
#include <cstdio>
#include <cstring>

#include <boost/scoped_array.hpp>

namespace ecell4
{

namespace gillespie
{

void GillespieSimulator::increment_molecules(const Species& sp)
{
    world_->add_molecules(sp, 1);

    for (boost::ptr_vector<ReactionRuleEvent>::iterator i(events_.begin());
        i != events_.end(); ++i)
    {
        (*i).inc(sp);
    }
}


void GillespieSimulator::decrement_molecules(const Species& sp)
{
    world_->remove_molecules(sp, 1);

    for (boost::ptr_vector<ReactionRuleEvent>::iterator i(events_.begin());
        i != events_.end(); ++i)
    {
        (*i).dec(sp);
    }
}

bool GillespieSimulator::__draw_next_reaction(void)
{
    std::vector<double> a(events_.size());
    const Real V(world_->volume());
    for (unsigned int idx(0); idx < events_.size(); ++idx)
    {
        // events_[idx].initialize(world_.get());
        a[idx] = events_[idx].propensity();
    }

    const double atot(std::accumulate(a.begin(), a.end(), double(0.0)));
    if (atot == 0.0)
    {
        // Any reactions cannot occur.
        this->dt_ = inf;
        return true;
    }

    const double rnd1(rng()->uniform(0, 1));
    const double dt(gsl_sf_log(1.0 / rnd1) / double(atot));
    const double rnd2(rng()->uniform(0, atot));

    int u(-1);
    double acc(0.0);
    const int len_a(a.size());
    do
    {
        u++;
        acc += a[u];
    } while (acc < rnd2 && u < len_a - 1);

    if (len_a == u)
    {
        // Any reactions cannot occur.
        this->dt_ = inf;
        return true;
    }

    next_reaction_rule_ = events_[u].reaction_rule();
    next_reaction_ = events_[u].draw();
    if (next_reaction_.k() <= 0.0)
    {
        this->dt_ += dt; // skip a reaction
        return false;
    }

    this->dt_ += dt;
    return true;
}

void GillespieSimulator::draw_next_reaction(void)
{
    if (events_.size() == 0)
    {
        this->dt_ = inf;
        return;
    }

    this->dt_ = 0.0;

    while (!__draw_next_reaction())
    {
        ; // pass
    }
}

void GillespieSimulator::step(void)
{
    last_reactions_.clear();

    if (this->dt_ == inf)
    {
        // Any reactions cannot occur.
        return;
    }

    const Real t0(t()), dt0(dt());

    if (dt0 == 0.0 || next_reaction_.k() <= 0.0)
    {
        // Any reactions cannot occur.
        return;
    }

    // Reaction[u] occurs.
    for (ReactionRule::reactant_container_type::const_iterator
        it(next_reaction_.reactants().begin());
        it != next_reaction_.reactants().end(); ++it)
    {
        decrement_molecules(*it);
    }

    for (ReactionRule::product_container_type::const_iterator
        it(next_reaction_.products().begin());
        it != next_reaction_.products().end(); ++it)
    {
        increment_molecules(*it);
    }

    this->set_t(t0 + dt0);
    num_steps_++;

    last_reactions_.push_back(std::make_pair(next_reaction_rule_, reaction_info_type(t(), next_reaction_.reactants(), next_reaction_.products())));

    this->draw_next_reaction();
}

bool GillespieSimulator::step(const Real &upto)
{
    if (upto <= t())
    {
        return false;
    }

    if (upto >= next_time())
    {
        step();
        return true;
    }
    else
    {
        // no reaction occurs
        // set_dt(next_time() - upto);
        set_t(upto);
        last_reactions_.clear();
        draw_next_reaction();
        return false;
    }
}

void GillespieSimulator::initialize(void)
{
    const Model::reaction_rule_container_type&
        reaction_rules(model_->reaction_rules());

    events_.clear();
    for (Model::reaction_rule_container_type::const_iterator
        i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);

        if (rr.reactants().size() == 0)
        {
            events_.push_back(new ZerothOrderReactionRuleEvent(this, rr));
        }
        else if (rr.reactants().size() == 1)
        {
            events_.push_back(new FirstOrderReactionRuleEvent(this, rr));
        }
        else if (rr.reactants().size() == 2)
        {
            events_.push_back(new SecondOrderReactionRuleEvent(this, rr));
        }
        else
        {
            throw NotSupported("not supported yet.");
        }

        events_.back().initialize();
    }

    this->draw_next_reaction();
}

Real GillespieSimulator::dt(void) const
{
    return this->dt_;
}

} // gillespie

} // ecell4
