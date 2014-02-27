#include "LatticeSimulator.hpp"


namespace ecell4
{

namespace lattice
{

void LatticeSimulator::initialize()
{
    //(*world_).initialize();
    dt_ = 0.1;
    // dt = (4 * R^2) / (6 * D) for each species

    std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        const Species species(*itr);
        if (boost::lexical_cast<Real>(species.get_attribute("D")) == 0.)
            continue;

        const boost::shared_ptr<EventScheduler::Event> event(create_event(species));
        scheduler_.add(event);
    }

    is_initialized_ = true;
}

boost::shared_ptr<EventScheduler::Event> LatticeSimulator::create_event(
        const Species& species)
{
    boost::shared_ptr<EventScheduler::Event> event(new Event(this, species));
    return event;
}

void LatticeSimulator::attempt_reaction_(Coord from_coord, Coord to_coord)
{
    // TODO
    const Species from_species(world_->get_molecular_type(from_coord)->species()),
                  to_species(world_->get_molecular_type(to_coord)->species());
    const std::vector<ReactionRule> rules(model_->query_reaction_rules(
                from_species, to_species));
    if (rules.size() != 1)
    {
        return;
    }
    const ReactionRule rule(rules.at(0));
    const Real k(rule.k());
    const ReactionRule::reactant_container_type& reactants(rule.reactants());
    if (reactants.count(from_species) != 1 ||
            reactants.count(to_species) != 1 ||
            reactants.size() != 2)
    {
        return;
    }
    const ReactionRule::product_container_type& products(rule.products());
    if (products.size() != 1)
    {
        return;
    }
    Species product_species(*(products.begin()));
    world_->remove_molecule(from_coord);
    world_->remove_molecule(to_coord);
    world_->add_molecule(product_species, to_coord);
}

void LatticeSimulator::step()
{
    if (!is_initialized_)
    {
        initialize();
    }

    EventScheduler::value_type const& top(scheduler_.top());
    Real time(top.second->time());
    top.second->fire(); // top.second->time_ is updated in fire()
    world_->set_t(time);
    scheduler_.update(top);
}

bool LatticeSimulator::step(const Real& upto)
{
    if (!is_initialized_)
    {
        initialize();
    }

    /*
    Integer count((Integer)(upto / dt()) - num_steps());
    for (Integer i(0); i < count; ++i)
    {
        step();
    }
    */

    if (upto < t())
    {
        return false;
    }

    EventScheduler::value_type next_event;
    while ((next_event = scheduler_.top()).second->time() < upto)
    {
        Real time(next_event.second->time());
        next_event.second->fire(); // top.second->time_ is updated in fire()
        world_->set_t(time);
        scheduler_.update(next_event);
    }

    return true;
}

} // lattice

} // ecell4
