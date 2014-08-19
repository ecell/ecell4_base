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

Integer GillespieSimulator::num_molecules(const Species& sp)
{
    SpeciesExpressionMatcher sexp(sp);
    Integer num_tot(0);
    const std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator i(species.begin());
        i != species.end(); ++i)
    {
        const Integer num(sexp.count(*i));
        if (num > 0)
        {
            num_tot += num * world_->num_molecules(*i); //XXX: num_molecules_exact
        }
    }
    return num_tot;
}

Integer GillespieSimulator::num_molecules(const Species& sp1, const Species& sp2)
{
    SpeciesExpressionMatcher sexp1(sp1), sexp2(sp2);
    Integer num_tot1(0), num_tot2(0), num_tot12(0);
    const std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator i(species.begin());
        i != species.end(); ++i)
    {
        const Integer num1(sexp1.count(*i));
        const Integer num2(sexp2.count(*i));

        if (num1 > 0 || num2 > 0)
        {
            const Integer num(world_->num_molecules(*i)); //XXX: num_molecules_exact
            const Integer tmp(num1 * num);
            num_tot1 += tmp;
            num_tot2 += num2 * num;
            num_tot12 += num2 * tmp;
        }
    }
    return num_tot1 * num_tot2 - num_tot12;
}

std::pair<ReactionRule::reactant_container_type, Integer>
GillespieSimulator::draw_exact_reactants(const Species& sp)
{
    const std::vector<Species> species(world_->list_species());

    SpeciesExpressionMatcher sexp(sp);
    Integer num_tot(0);
    std::vector<Real> tmp;
    for (std::vector<Species>::const_iterator j(species.begin());
        j != species.end(); ++j)
    {
        const Integer num(sexp.count(*j));
        if (num > 0)
        {
            num_tot += num * world_->num_molecules(*j); //XXX: num_molecules_exact
        }
        tmp.push_back(num_tot);
    }

    const Real rnd1(rng()->uniform(0.0, num_tot));
    std::vector<Real>::iterator
        itr(std::lower_bound(tmp.begin(), tmp.end(), rnd1));
    const Species& tgt(species[std::distance(tmp.begin(), itr)]);

    if (world_->num_molecules(tgt) == 0)
    {
        throw IllegalState("the number of reactant molecules must be non-zero.");
    }

    return std::make_pair(ReactionRule::reactant_container_type(1, tgt), sexp.count(tgt));
}

std::pair<ReactionRule::reactant_container_type, Integer>
GillespieSimulator::draw_exact_reactants(const Species& sp1, const Species& sp2)
{
    const std::vector<Species> species(world_->list_species());

    SpeciesExpressionMatcher sexp1(sp1), sexp2(sp2);

    Integer num_tot(0), cmb(0);

    std::vector<Real> tmp;
    for (std::vector<Species>::const_iterator j(species.begin());
        j != species.end(); ++j)
    {
        const Integer num1(sexp1.count(*j));
        if (num1 > 0)
        {
            num_tot += num1 * world_->num_molecules(*j); //XXX: num_molecules_exact
        }
        tmp.push_back(num_tot);
    }

    const Real rnd1(rng()->uniform(0.0, num_tot));
    std::vector<Real>::iterator
        itr1(std::lower_bound(tmp.begin(), tmp.end(), rnd1));
    const std::vector<Species>::difference_type
        idx1(std::distance(tmp.begin(), itr1));
    const Species& tgt1(species[idx1]);
    cmb = sexp1.count(tgt1);

    if (world_->num_molecules(tgt1) == 0)
    {
        throw IllegalState("the number of reactant molecules must be non-zero.");
    }

    tmp.clear();
    num_tot = 0;
    for (std::vector<Species>::const_iterator j(species.begin());
        j != species.end(); ++j)
    {
        const Integer num2(sexp2.count(*j));
        if (num2 > 0)
        {
            const Integer num(world_->num_molecules(*j)); //XXX: num_molecules_exact
            if (std::distance(species.begin(), j) != idx1)
            {
                num_tot += num2 * num;
            }
            else
            {
                num_tot += num2 * (num - 1);
            }
        }
        tmp.push_back(num_tot);
    }

    const Real rnd2(rng()->uniform(0.0, num_tot));
    std::vector<Real>::iterator
        itr2(std::lower_bound(tmp.begin(), tmp.end(), rnd2));
    const Species& tgt2(species[std::distance(tmp.begin(), itr2)]);
    cmb *= sexp2.count(tgt2);

    ReactionRule::reactant_container_type reactants(2);
    reactants[0] = tgt1;
    reactants[1] = tgt2;
    return std::make_pair(reactants, cmb);
}

ReactionRule GillespieSimulator::draw_exact_reaction(const ReactionRule& rr)
{
    ReactionRule::reactant_container_type reactants;
    Integer cmb;
    if (rr.reactants().size() == 1)
    {
        std::pair<ReactionRule::reactant_container_type, Integer> retval(draw_exact_reactants(rr.reactants()[0]));
        reactants = retval.first;
        cmb = retval.second;
    }
    else if (rr.reactants().size() == 2)
    {
        std::pair<ReactionRule::reactant_container_type, Integer> retval(draw_exact_reactants(rr.reactants()[0], rr.reactants()[1]));
        reactants = retval.first;
        cmb = retval.second;
    }
    else
    {
        throw NotSupported("not supported yet.");
    }

    std::vector<std::vector<Species> > possible_products(rrgenerate(rr, reactants));
    if (cmb <= 0 || cmb < possible_products.size())
    {
        return ReactionRule(); //XXX: error?
    }
    else if (possible_products.size() == 0)
    {
        return ReactionRule();
    }
    else if (cmb == 1)
    {
        // assert(possible_products.size() == 1);
        return ReactionRule(
            reactants, possible_products[0], rr.k());
    }
    else
    {
        const Integer rnd2(rng()->uniform_int(0, cmb - 1));
        if (rnd2 >= possible_products.size())
        {
            return ReactionRule();
        }

        return ReactionRule(
            reactants, possible_products[rnd2], rr.k());
    }
}

Real GillespieSimulator::calculate_propensity(const ReactionRule& rr)
{
    if (rr.reactants().size() == 1)
    {
        return rr.k() * num_molecules(rr.reactants()[0]);
    }
    else if (rr.reactants().size() == 2)
    {
        const Real V(world_->volume());
        return (rr.k() * num_molecules(rr.reactants()[0], rr.reactants()[1])
            / world_->volume());
    }
    else
    {
        throw NotSupported("not supported yet.");
    }
}

bool GillespieSimulator::__draw_next_reaction(void)
{
    const Model::reaction_rule_container_type&
        reaction_rules(model_->reaction_rules());

    std::vector<double> a(reaction_rules.size());
    for (unsigned int idx(0); idx < reaction_rules.size(); ++idx)
    {
        a[idx] = calculate_propensity(reaction_rules[idx]);
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

    next_reaction_ = draw_exact_reaction(reaction_rules[u]);
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
    if (model_->reaction_rules().size() == 0)
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
        world_->remove_molecules(*it, 1);
    }

    for (ReactionRule::product_container_type::const_iterator
        it(next_reaction_.products().begin());
        it != next_reaction_.products().end(); ++it)
    {
        world_->add_molecules(*it, 1);
    }

    last_reactions_.clear();
    last_reactions_.push_back(next_reaction_);

    this->set_t(t0 + dt0);
    num_steps_++;

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
        // last_reactions_.clear();
        draw_next_reaction();
        return false;
    }
}

void GillespieSimulator::initialize(void)
{
    this->draw_next_reaction();
}

void GillespieSimulator::set_t(const Real &t)
{
    this->world_->set_t(t);
}

Real GillespieSimulator::t(void) const
{
    return this->world_->t();
}

Real GillespieSimulator::dt(void) const
{
    return this->dt_;
}

std::vector<ReactionRule> GillespieSimulator::last_reactions() const
{
    return last_reactions_;
}

} // gillespie

} // ecell4
