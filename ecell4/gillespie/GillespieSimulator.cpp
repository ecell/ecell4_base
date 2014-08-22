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

void GillespieSimulator::add_molecules(const Species& sp)
{
    world_->add_molecules(sp, 1);

    for (boost::ptr_vector<ReactionRuleEvent>::iterator i(events_.begin());
        i != events_.end(); ++i)
    {
        (*i).inc(sp);
    }
}


void GillespieSimulator::remove_molecules(const Species& sp)
{
    world_->remove_molecules(sp, 1);

    for (boost::ptr_vector<ReactionRuleEvent>::iterator i(events_.begin());
        i != events_.end(); ++i)
    {
        (*i).dec(sp);
    }
}

// void GillespieSimulator::calculate_stoichiometries()
// {
//     const Model::reaction_rule_container_type&
//         reaction_rules(model_->reaction_rules());
// 
//     stoichiometries_.clear();
//     for (Model::reaction_rule_container_type::const_iterator
//         i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
//     {
//         const ReactionRule::reactant_container_type& reactants((*i).reactants());
//         if (reactants.size() == 1)
//         {
//             stoichiometries_.push_back(get_stoichiometry(reactants[0]));
//         }
//         else if (reactants.size() == 2)
//         {
//             stoichiometries_.push_back(get_stoichiometry(reactants[0], reactants[1]));
//         }
//         else
//         {
//             throw NotSupported("not supported yet.");
//         }
//     }
// }
// 
// void GillespieSimulator::append_stoichiometries(const Species& sp)
// {
//     const Model::reaction_rule_container_type&
//         reaction_rules(model_->reaction_rules());
// 
//     for (unsigned int i(0); i < reaction_rules.size(); ++i)
//     {
//         const ReactionRule::reactant_container_type&
//             reactants(reaction_rules[i].reactants());
//         if (reactants.size() == 1)
//         {
//             const Integer coef(SpeciesExpressionMatcher(reactants[0]).count(sp));
//             if (coef > 0)
//             {
//                 stoichiometries_[i].push_back(stoichiometry(sp, coef));
//             }
//         }
//         else if (reactants.size() == 2)
//         {
//             const Integer coef1(SpeciesExpressionMatcher(reactants[0]).count(sp));
//             const Integer coef2(SpeciesExpressionMatcher(reactants[1]).count(sp));
//             if (coef1 > 0 || coef2 > 0)
//             {
//                 stoichiometries_[i].push_back(stoichiometry(sp, coef1, coef2));
//             }
//         }
//         else
//         {
//             throw NotSupported("not supported yet.");
//         }
//     }
// }
// 
// GillespieSimulator::stoichiometry_container_type
// GillespieSimulator::get_stoichiometry(const Species& sp)
// {
//     const std::vector<Species> species(world_->list_species());
//     stoichiometry_container_type retval;
// 
//     SpeciesExpressionMatcher sexp(sp);
//     for (std::vector<Species>::const_iterator i(species.begin());
//         i != species.end(); ++i)
//     {
//         const Integer coef(sexp.count(*i));
//         if (coef > 0)
//         {
//             retval.push_back(stoichiometry(*i, coef));
//         }
//     }
//     return retval;
// }
// 
// GillespieSimulator::stoichiometry_container_type
// GillespieSimulator::get_stoichiometry(const Species& sp1, const Species& sp2)
// {
//     const std::vector<Species> species(world_->list_species());
//     stoichiometry_container_type retval;
// 
//     SpeciesExpressionMatcher sexp1(sp1), sexp2(sp2);
//     for (std::vector<Species>::const_iterator i(species.begin());
//         i != species.end(); ++i)
//     {
//         const Integer coef1(sexp1.count(*i));
//         const Integer coef2(sexp2.count(*i));
//         if (coef1 > 0 || coef2 > 0)
//         {
//             retval.push_back(stoichiometry(*i, coef1, coef2));
//         }
//     }
//     return retval;
// }
// 
// Integer GillespieSimulator::num_molecules(
//     const Model::reaction_rule_container_type::size_type& u)
// {
//     const stoichiometry_container_type& retval(stoichiometries_[u]);
//     const ReactionRule::reactant_container_type& reactants(
//         model_->reaction_rules()[u].reactants());
// 
//     if (reactants.size() == 1)
//     {
//         Integer num_tot(0);
//         for (stoichiometry_container_type::const_iterator
//             i(retval.begin()); i != retval.end(); ++i)
//         {
//             num_tot += (*i).coef1 * world_->num_molecules_exact((*i).species);
//         }
//         return num_tot;
//     }
//     else if (reactants.size() == 2)
//     {
//         Integer num_tot1(0), num_tot2(0), num_tot12(0);
//         for (stoichiometry_container_type::const_iterator
//             i(retval.begin()); i != retval.end(); ++i)
//         {
//             const Integer num(world_->num_molecules_exact((*i).species));
//             const Integer tmp(num * (*i).coef1);
//             num_tot1 += tmp;
//             num_tot2 += num * (*i).coef2;
//             num_tot12 += tmp * (*i).coef2;
//         }
//         return num_tot1 * num_tot2 - num_tot12;
//     }
//     else
//     {
//         throw NotSupported("not supported yet.");
//     }
// }
// 
// std::pair<ReactionRule::reactant_container_type, Integer>
// GillespieSimulator::draw_exact_reactants(
//     const Model::reaction_rule_container_type::size_type& u)
// {
//     ReactionRule::reactant_container_type reactants(
//         model_->reaction_rules()[u].reactants());
//     if (reactants.size() == 1)
//     {
//         return draw_exact_reactants(reactants[0], stoichiometries_[u]);
//     }
//     else if (reactants.size() == 2)
//     {
//         return draw_exact_reactants(
//             reactants[0], reactants[1], stoichiometries_[u]);
//     }
//     else
//     {
//         throw NotSupported("not supported yet.");
//     }
// }
// 
// std::pair<ReactionRule::reactant_container_type, Integer>
// GillespieSimulator::draw_exact_reactants(
//     const Species& sp, const stoichiometry_container_type& retval)
// {
//     if (retval.size() == 1)
//     {
//         return std::make_pair(
//             ReactionRule::reactant_container_type(1, retval[0].species),
//             retval[0].coef1);
//     }
// 
//     Integer num_tot(0);
//     std::vector<Real> cum;
//     cum.reserve(retval.size());
//     for (stoichiometry_container_type::const_iterator
//         i(retval.begin()); i != retval.end(); ++i)
//     {
//         num_tot += (*i).coef1 * world_->num_molecules_exact((*i).species);
//         cum.push_back(num_tot);
//     }
// 
//     assert(retval.size() > 0);
//     const Real rnd1(rng()->uniform(0.0, num_tot));
//     unsigned int idx(std::distance(
//         cum.begin(), std::lower_bound(cum.begin(), cum.end(), rnd1)));
//     assert(idx < retval.size());
//     const Species& tgt(retval[idx].species);
//     const Integer& coef(retval[idx].coef1);
// 
//     // assert(world_->num_molecules_exact(tgt) > 0);
//     return std::make_pair(ReactionRule::reactant_container_type(1, tgt), coef);
// }
// 
// std::pair<ReactionRule::reactant_container_type, Integer>
// GillespieSimulator::draw_exact_reactants(
//     const Species& sp1, const Species& sp2, const stoichiometry_container_type& retval)
// {
//     if (retval.size() == 2)
//     {
//         if (retval[0].coef1 == 0)
//         {
//             if (retval[1].coef2 == 0)
//             {
//                 ReactionRule::reactant_container_type reactants(2);
//                 reactants[0] = retval[1].species;
//                 reactants[1] = retval[0].species;
//                 return std::make_pair(reactants,
//                     retval[0].coef2 * retval[1].coef1);
//             }
//         }
//         else if (retval[1].coef1 == 0)
//         {
//             if (retval[0].coef2 == 0)
//             {
//                 ReactionRule::reactant_container_type reactants(2);
//                 reactants[0] = retval[0].species;
//                 reactants[1] = retval[1].species;
//                 return std::make_pair(reactants,
//                     retval[0].coef1 * retval[1].coef2);
//             }
//         }
//     }
// 
//     Integer num_tot1(0), num_tot2(0);
//     std::vector<Real> cum1, cum2;
//     cum1.reserve(retval.size());
//     cum2.reserve(retval.size());
//     for (stoichiometry_container_type::const_iterator
//         i(retval.begin()); i != retval.end(); ++i)
//     {
//         const Integer num(world_->num_molecules_exact((*i).species));
//         num_tot1 += (*i).coef1 * num;
//         num_tot2 += (*i).coef2 * num;
//         cum1.push_back(num_tot1);
//         cum2.push_back(num_tot2);
//     }
// 
//     assert(retval.size() > 0);
//     const Real rnd1(rng()->uniform(0.0, num_tot1));
//     unsigned int idx1(std::distance(
//         cum1.begin(), std::lower_bound(cum1.begin(), cum1.end(), rnd1)));
//     assert(idx1 < retval.size());
//     const Species& tgt1(retval[idx1].species);
//     const Integer& coef1(retval[idx1].coef1);
//     // assert(world_->num_molecules_exact(tgt1) > 0);
// 
//     const Integer& coef12(retval[idx1].coef2);
//     if (coef12 > 0)
//     {
//         num_tot2 -= coef12;
//         std::transform(cum2.begin() + idx1, cum2.end(),
//             cum2.begin() + idx1, std::bind2nd(std::minus<Real>(), coef12));
//     }
// 
//     const Real rnd2(rng()->uniform(0.0, num_tot2));
//     unsigned int idx2(std::distance(
//         cum2.begin(), std::lower_bound(cum2.begin(), cum2.end(), rnd2)));
//     assert(idx2 < retval.size());
//     const Species& tgt2(retval[idx2].species);
//     const Integer& coef2(retval[idx2].coef2);
//     // assert(world_->num_molecules_exact(tgt2) > 0);
// 
//     ReactionRule::reactant_container_type reactants(2);
//     reactants[0] = tgt1;
//     reactants[1] = tgt2;
//     return std::make_pair(reactants, coef1 * coef2);
// }
// 
// ReactionRule GillespieSimulator::draw_exact_reaction(
//     const Model::reaction_rule_container_type::size_type& u)
// {
//     const ReactionRule& rr(model_->reaction_rules()[u]);
//     std::pair<ReactionRule::reactant_container_type, Integer>
//         retval(draw_exact_reactants(u));
//     ReactionRule::reactant_container_type& reactants(retval.first);
//     const Integer cmb(retval.second);
// 
//     std::vector<std::vector<Species> > possible_products(rrgenerate(rr, reactants));
//     if (cmb <= 0 || cmb < possible_products.size())
//     {
//         return ReactionRule(); //XXX: error?
//     }
//     else if (possible_products.size() == 0)
//     {
//         return ReactionRule();
//     }
//     else if (cmb == 1)
//     {
//         // assert(possible_products.size() == 1);
//         return ReactionRule(
//             reactants, possible_products[0], rr.k());
//     }
//     else
//     {
//         const Integer rnd2(rng()->uniform_int(0, cmb - 1));
//         if (rnd2 >= possible_products.size())
//         {
//             return ReactionRule();
//         }
// 
//         return ReactionRule(
//             reactants, possible_products[rnd2], rr.k());
//     }
// }

// Real GillespieSimulator::calculate_propensity(
//     const Model::reaction_rule_container_type::size_type& u)
// {
//     const ReactionRule& rr(model_->reaction_rules()[u]);
//     const Real V(world_->volume());
//     Real prop(rr.k() * num_molecules(u));
//     for (unsigned int i(0); i < rr.reactants().size() - 1; ++i)
//     {
//         prop /= V;
//     }
//     return prop;
// }

bool GillespieSimulator::__draw_next_reaction(void)
{
    std::vector<double> a(events_.size());
    const Real V(world_->volume());
    for (unsigned int idx(0); idx < events_.size(); ++idx)
    {
        // events_[idx].initialize(world_.get());
        a[idx] = events_[idx].propensity(V);
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

    next_reaction_ = events_[u].draw(world_.get());
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
        remove_molecules(*it);
    }

    for (ReactionRule::product_container_type::const_iterator
        it(next_reaction_.products().begin());
        it != next_reaction_.products().end(); ++it)
    {
        add_molecules(*it);
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
    const Model::reaction_rule_container_type&
        reaction_rules(model_->reaction_rules());

    events_.clear();
    for (Model::reaction_rule_container_type::const_iterator
        i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);

        if (rr.reactants().size() == 1)
        {
            events_.push_back(new FirstOrderReactionRuleEvent(rr));
        }
        else if (rr.reactants().size() == 2)
        {
            events_.push_back(new SecondOrderReactionRuleEvent(rr));
        }
        else
        {
            throw NotSupported("not supported yet.");
        }

        events_.back().initialize(world_.get());
    }

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
