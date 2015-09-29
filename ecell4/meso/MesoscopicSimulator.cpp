#include <numeric>
#include <vector>
#include <gsl/gsl_sf_log.h>

#include <cstring>
#include <sstream>
#include <cstdio>
#include <cstring>

#include <boost/scoped_array.hpp>

#include "MesoscopicSimulator.hpp"


namespace ecell4
{

namespace meso
{

void MesoscopicSimulator::increment_molecules(const Species& sp, const coordinate_type& c)
{
    world_->add_molecules(sp, 1, c);

    for (boost::ptr_vector<ReactionRuleProxyBase>::iterator i(proxies_.begin());
        i != proxies_.end(); ++i)
    {
        (*i).inc(sp, c);
    }
}

void MesoscopicSimulator::decrement_molecules(const Species& sp, const coordinate_type& c)
{
    world_->remove_molecules(sp, 1, c);

    for (boost::ptr_vector<ReactionRuleProxyBase>::iterator i(proxies_.begin());
        i != proxies_.end(); ++i)
    {
        (*i).dec(sp, c);
    }
}

std::pair<Real, std::pair<ReactionRule, MesoscopicSimulator::coordinate_type> >
MesoscopicSimulator::draw_next_reaction(const coordinate_type& c)
{
    std::vector<double> a(proxies_.size());
    for (unsigned int idx(0); idx < proxies_.size(); ++idx)
    {
        a[idx] = proxies_[idx].propensity(c);
    }

    const double atot(std::accumulate(a.begin(), a.end(), double(0.0)));
    if (atot == 0.0)
    {
        // Any reactions cannot occur.
        return std::make_pair(inf, std::make_pair(ReactionRule(), c));
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
        return std::make_pair(inf, std::make_pair(ReactionRule(), c));
    }

    return std::make_pair(dt, proxies_[u].draw(c));
}

void MesoscopicSimulator::interrupt_all(const Real& t)
{
    EventScheduler::events_range events(scheduler_.events());
    for (EventScheduler::events_range::iterator itr(events.begin());
            itr != events.end(); ++itr)
    {
        (*itr).second->interrupt(t);
        scheduler_.update(*itr);
    }
}

void MesoscopicSimulator::step(void)
{
    if (this->dt() == inf)
    {
        // Any reactions cannot occur.
        return;
    }

    interrupted_ = event_ids_.size();
    EventScheduler::value_type const& top(scheduler_.top());
    const Real tnext(top.second->time());
    top.second->fire(); // top.second->time_ is updated in fire()
    this->set_t(tnext);
    scheduler_.update(top);

    if (interrupted_ < event_ids_.size())
    {
        EventScheduler::identifier_type evid(event_ids_[interrupted_]);
        boost::shared_ptr<EventScheduler::Event> ev(scheduler_.get(evid));
        ev->interrupt(t());
        scheduler_.update(std::make_pair(evid, ev));
    }

    // EventScheduler::value_type top(scheduler_.pop());
    // const Real tnext(top.second->time());
    // top.second->fire(); // top.second->time_ is updated in fire()
    // this->set_t(tnext);
    // // EventScheduler::events_range events(scheduler_.events());
    // // for (EventScheduler::events_range::iterator itr(events.begin());
    // //         itr != events.end(); ++itr)
    // // {
    // //     (*itr).second->interrupt(tnext);
    // //     scheduler_.update(*itr);
    // // }
    // scheduler_.add(top.second);

    num_steps_++;
}

bool MesoscopicSimulator::step(const Real &upto)
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
        // nothing happens
        // set_dt(next_time() - upto);
        set_t(upto);
        last_reactions_.clear();
        // interrupt_all(upto);  //XXX: Is this really needed?
        return false;
    }
}

void MesoscopicSimulator::initialize(void)
{
    const Model::reaction_rule_container_type&
        reaction_rules(model_->reaction_rules());

    proxies_.clear();
    for (Model::reaction_rule_container_type::const_iterator
        i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);

        if (rr.reactants().size() == 0)
        {
            proxies_.push_back(new ZerothOrderReactionRuleProxy(this, rr));
        }
        else if (rr.reactants().size() == 1)
        {
            proxies_.push_back(new FirstOrderReactionRuleProxy(this, rr));
        }
        else if (rr.reactants().size() == 2)
        {
            if (world_->has_structure(rr.reactants()[0]))
            {
                proxies_.push_back(new StructureSecondOrderReactionRuleProxy(this, rr, 0));
            }
            else if (world_->has_structure(rr.reactants()[1]))
            {
                proxies_.push_back(new StructureSecondOrderReactionRuleProxy(this, rr, 1));
            }
            else
            {
                proxies_.push_back(new SecondOrderReactionRuleProxy(this, rr));
            }
        }
        else
        {
            throw NotSupported("not supported yet.");
        }

        proxies_.back().initialize();
    }

    const std::vector<Species>& species(model_->species_attributes());
    // const std::vector<Species>& species(world_->species());
    for (std::vector<Species>::const_iterator i(species.begin());
        i != species.end(); ++i)
    {
        proxies_.push_back(new DiffusionProxy(this, *i));
        proxies_.back().initialize();
    }

    scheduler_.clear();
    event_ids_.resize(world_->num_subvolumes());
    for (Integer i(0); i < world_->num_subvolumes(); ++i)
    {
        event_ids_[i] =
            scheduler_.add(boost::shared_ptr<EventScheduler::Event>(
                new SubvolumeEvent(this, i, t())));
    }
}

Real MesoscopicSimulator::dt(void) const
{
    return next_time() - t();
}

Real MesoscopicSimulator::next_time(void) const
{
    return scheduler_.next_time();
}

std::vector<ReactionRule> MesoscopicSimulator::last_reactions() const
{
    return last_reactions_;
}

void MesoscopicSimulator::set_last_reaction(const ReactionRule& rr)
{
    last_reactions_.clear();
    last_reactions_.push_back(rr);
}

} // meso

} // ecell4
