#include "GillespieSimulator.hpp"
#include <numeric>
#include <vector>
#include <gsl/gsl_sf_log.h>

namespace ecell4
{

namespace gillespie
{

void GillespieSimulator::calc_next_reaction_(void)
{
    // reset
    this->dt_ = inf;

    const NetworkModel::reaction_rule_container_type &possible_reaction_rules =
        this->model_->reaction_rules();
    if (possible_reaction_rules.size() == 0)
    {
        this->dt_ = inf;
        return;
    }
    std::vector<double> a(possible_reaction_rules.size());
    for (unsigned int idx(0); idx < possible_reaction_rules.size(); idx++)
    {
        a[idx] = possible_reaction_rules[idx].k() * this->world_->volume();
        const ReactionRule::reactant_container_type &reactants =
            possible_reaction_rules[idx].reactants();
        for (ReactionRule::reactant_container_type::iterator
                 it = reactants.begin();
             it != reactants.end(); it++)
        {
            a[idx] *= this->world_->num_molecules(*it) / this->world_->volume();
        }
    }
    double a_total(std::accumulate(a.begin(), a.end(), double(0.0)));
    if(a_total == 0.0)
    {
        // Any reactions cannot occur.
        this->dt_ = inf;
        return;
    }

    double rnd_num1(this->rng()->uniform(0, 1));
    double dt(gsl_sf_log(1.0 / rnd_num1) / double(a_total));
    double rnd_num2(this->rng()->uniform(0, 1) * a_total);

    int u(-1);
    double acc(0.0);
    int len(a.size());
    do
    {
        u++;
        acc += a[u];
    } while (acc < rnd_num2 && u < len -1);

    if (len == u)
    {
        // Any reactions cannot occur.
        this->dt_ = inf;
        return;
    }

    // save.
    this->next_reaction_num_ = u;
    this->dt_ = dt;
}

void GillespieSimulator::step(void)
{
    if (this->dt_ == inf)
    {
        // Any reactions cannot occur.
        return;
    }
    const NetworkModel::reaction_rule_container_type &possible_reaction_rules =
        this->model_->reaction_rules();

    int u = this->next_reaction_num_;
    Real dt = this->dt_;

    if (dt == 0.0 || u < 0)
    {
        // Any reactions cannot occur.
        return;
    }

    //Reaction[u] occurs.
    for (ReactionRule::reactant_container_type::iterator
            it(possible_reaction_rules[u].reactants().begin());
        it != possible_reaction_rules[u].reactants().end();
        it++)
    {
        int one(1);
        this->world_->remove_molecules(*it, one);
    }
    for (ReactionRule::product_container_type::iterator
             it(possible_reaction_rules[u].products().begin());
        it != possible_reaction_rules[u].products().end();
        it++)
    {
        int one(1);
        this->world_->add_molecules(*it, one);
    }
    this->world_->set_t(this->world_->t() + dt);
    this->num_steps_++;

    this->calc_next_reaction_();
}

bool GillespieSimulator::step(Real const &upto)
{
    Real const t0(t()), dt0(dt());
    Real const next_time(t0 + dt0);

    if (upto <= t0)
    {
        return false;
    }

    if (upto >= next_time)
    {
        this->step();
        return true;
    }
    else
    {
        // no reaction occurs
        this->set_t(upto);
        this->initialize();
        return false;
    }
}

void GillespieSimulator::initialize(void)
{
    this->calc_next_reaction_();
}

void GillespieSimulator::set_t(Real const &t)
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

Integer GillespieSimulator::num_steps(void) const
{
    return this->num_steps_;
}

} // gillespie

} // ecell4
