#include "BDSimulator.hpp"

#include <cstring>

namespace ecell4
{

namespace bd
{

void BDSimulator::attempt_synthetic_reaction(const ReactionRule& rr)
{
    assert(rr.reactants().size() == 0);

    const Real pacc = rr.k() * (*world_).volume() * dt();  //XXX: dt must be small enough
    assert(pacc <= 1);
    // const Real pacc = 1.0 - exp(-rr.k() * (*world_).volume() * dt());
    const Real rnd = rng()->uniform(0, 1);
    if (rnd >= pacc)
    {
        return;
    }

    reaction_info_type::container_type new_particle_ids;

    {
        for (ReactionRule::product_container_type::const_iterator i(rr.products().begin());
            i != rr.products().end(); ++i)
        {
            const Real3 newpos(
                rng()->uniform(0, (*world_).edge_lengths()[0]),
                rng()->uniform(0, (*world_).edge_lengths()[1]),
                rng()->uniform(0, (*world_).edge_lengths()[2]));
            std::pair<std::pair<ParticleID, Particle>, bool>
                ret = (*world_).new_particle((*i), newpos);
            if (ret.second)
            {
                new_particle_ids.push_back(ret.first);
            }
            else
            {
                for (reaction_info_type::container_type::const_iterator
                    j(new_particle_ids.begin()); j != new_particle_ids.end(); ++j)
                {
                    (*world_).remove_particle((*j).first);
                }
                return;
            }
        }
    }

    reaction_info_type ri(t() + dt(), reaction_info_type::container_type(), new_particle_ids);
    last_reactions_.push_back(std::make_pair(rr, ri));
}

void BDSimulator::step()
{
    last_reactions_.clear();

    for (Model::reaction_rule_container_type::const_iterator i((*model_).reaction_rules().begin());
        i != (*model_).reaction_rules().end(); ++i)
    {
        if ((*i).reactants().size() == 0)
        {
            attempt_synthetic_reaction(*i);
        }
    }

    {
        BDPropagator propagator(*model_, *world_, *rng(), dt(), last_reactions_);
        while (propagator())
        {
            ; // do nothing here
        }
    }

    set_t(t() + dt());
    num_steps_++;
}

bool BDSimulator::step(const Real& upto)
{
    const Real t0(t()), dt0(dt()), tnext(next_time());

    if (upto <= t0)
    {
        return false;
    }

    if (upto >= tnext)
    {
        step();
        return true;
    }
    else
    {
        dt_ = upto - t0;
        step();
        dt_ = dt0;
        return false;
    }
}

} // bd

} // ecell4
