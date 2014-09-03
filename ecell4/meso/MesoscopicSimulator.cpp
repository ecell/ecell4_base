#include "MesoscopicSimulator.hpp"
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

namespace meso
{

void MesoscopicSimulator::increment_molecules(const Species& sp, const coordinate_type& c)
{
    world_->add_molecules(sp, 1, c);

    for (boost::ptr_vector<ReactionRuleProxy>::iterator i(proxies_.begin());
        i != proxies_.end(); ++i)
    {
        (*i).inc(sp, c);
    }

    for (boost::ptr_vector<DiffusionProxy>::iterator i(diffusions_.begin());
        i != diffusions_.end(); ++i)
    {
        (*i).inc(sp, c);
    }
}

void MesoscopicSimulator::decrement_molecules(const Species& sp, const coordinate_type& c)
{
    world_->remove_molecules(sp, 1, c);

    for (boost::ptr_vector<ReactionRuleProxy>::iterator i(proxies_.begin());
        i != proxies_.end(); ++i)
    {
        (*i).dec(sp, c);
    }

    for (boost::ptr_vector<DiffusionProxy>::iterator i(diffusions_.begin());
        i != diffusions_.end(); ++i)
    {
        (*i).dec(sp, c);
    }
}

std::pair<Real, std::pair<Species, MesoscopicSimulator::coordinate_type> >
MesoscopicSimulator::draw_next_diffusion(const coordinate_type& c)
{
    std::vector<double> b(diffusions_.size());
    for (unsigned int idx(0); idx < diffusions_.size(); ++idx)
    {
        b[idx] = diffusions_[idx].propensity(c);
    }

    const double btot(std::accumulate(b.begin(), b.end(), double(0.0)));

    if (btot == 0.0)
    {
        // Any reactions cannot occur.
        return std::make_pair(inf, std::make_pair(Species(), c));
    }

    const double rnd1(rng()->uniform(0, 1));
    const double dt(gsl_sf_log(1.0 / rnd1) / double(btot));
    const double rnd2(rng()->uniform(0, btot));

    int u(-1);
    double acc(0.0);
    const int len_b(b.size());
    do
    {
        u++;
        acc += b[u];
    } while (acc < rnd2 && u < len_b - 1);

    if (len_b == u)
    {
        // Any reactions cannot occur.
        return std::make_pair(inf, std::make_pair(Species(), c));
    }

    const std::pair<Species, coordinate_type> retval(diffusions_[u].draw(c));
    assert(b[u] > 0);
    assert(world_->num_molecules_exact(retval.first, c) > 0);
    return std::make_pair(dt, retval);
}

std::pair<Real, ReactionRule> MesoscopicSimulator::draw_next_reaction(const coordinate_type& c)
{
    std::vector<double> a(proxies_.size());
    const Real V(world_->volume());
    for (unsigned int idx(0); idx < proxies_.size(); ++idx)
    {
        // proxies_[idx].initialize(world_.get());
        a[idx] = proxies_[idx].propensity(c);
    }

    const double atot(std::accumulate(a.begin(), a.end(), double(0.0)));
    if (atot == 0.0)
    {
        // Any reactions cannot occur.
        return std::make_pair(inf, ReactionRule());
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
        return std::make_pair(inf, ReactionRule());
    }

    const ReactionRule next_reaction(proxies_[u].draw(c));
    return std::make_pair(dt, next_reaction);
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
    const Real time(top.second->time());
    top.second->fire(); // top.second->time_ is updated in fire()
    this->set_t(time);
    scheduler_.update(top);

    if (interrupted_ < event_ids_.size())
    {
        EventScheduler::identifier_type evid(event_ids_[interrupted_]);
        boost::shared_ptr<EventScheduler::Event> ev(scheduler_.get(evid));
        ev->interrupt(t());
        scheduler_.update(std::make_pair(evid, ev));
    }

    // EventScheduler::value_type top(scheduler_.pop());
    // const Real time(top.second->time());
    // top.second->fire(); // top.second->time_ is updated in fire()
    // this->set_t(time);
    // // EventScheduler::events_range events(scheduler_.events());
    // // for (EventScheduler::events_range::iterator itr(events.begin());
    // //         itr != events.end(); ++itr)
    // // {
    // //     (*itr).second->interrupt(time);
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
        // no reaction occurs
        // set_dt(next_time() - upto);
        set_t(upto);
        // last_reactions_.clear();
        interrupt_all(upto);
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

        if (rr.reactants().size() == 1)
        {
            proxies_.push_back(new FirstOrderReactionRuleProxy(this, rr));
        }
        else if (rr.reactants().size() == 2)
        {
            proxies_.push_back(new SecondOrderReactionRuleProxy(this, rr));
        }
        else
        {
            throw NotSupported("not supported yet.");
        }

        proxies_.back().initialize();
    }

    const std::vector<Species>& species(world_->species());
    diffusions_.clear();
    for (std::vector<Species>::const_iterator i(species.begin());
        i != species.end(); ++i)
    {
        diffusions_.push_back(new DiffusionProxy(this, *i));
        diffusions_.back().initialize();
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

void MesoscopicSimulator::set_t(const Real &t)
{
    this->world_->set_t(t);
}

Real MesoscopicSimulator::t(void) const
{
    return this->world_->t();
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
