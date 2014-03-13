#include "LatticeSimulator.hpp"


namespace ecell4
{

namespace lattice
{

void LatticeSimulator::initialize()
{
    scheduler_.clear();
    std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        const Species species(*itr);
        const boost::shared_ptr<EventScheduler::Event> event(
                create_step_event(species));
        scheduler_.add(event);
    }

    NetworkModel::reaction_rule_container_type rules(model_->reaction_rules());
    for (NetworkModel::reaction_rule_container_type::iterator itr(rules.begin());
            itr != rules.end(); ++itr)
    {
        if ((*itr).reactants().size() == 1)
        {
            const boost::shared_ptr<EventScheduler::Event> event(
                    create_first_order_reaction_event(*itr));
            scheduler_.add(event);
        }
    }

    is_initialized_ = true;
}

boost::shared_ptr<EventScheduler::Event> LatticeSimulator::create_step_event(
        const Species& species)
{
    boost::shared_ptr<EventScheduler::Event> event(new StepEvent(this, species));
    return event;
}

boost::shared_ptr<EventScheduler::Event>
LatticeSimulator::create_first_order_reaction_event(const ReactionRule& reaction_rule)
{
    boost::shared_ptr<EventScheduler::Event> event(new FirstOrderReactionEvent(
                this, reaction_rule));
    return event;
}


void LatticeSimulator::attempt_reaction_(Coord from_coord, Coord to_coord)
{
    const Species from_species(world_->get_molecular_type(from_coord)->species()),
                  to_species(world_->get_molecular_type(to_coord)->species());
    const std::vector<ReactionRule> rules(model_->query_reaction_rules(
                from_species, to_species));

    const Real factor((from_species == to_species) ? 1 : 2);
    const Real Da(boost::lexical_cast<Real>(from_species.get_attribute("D"))),
               Db(boost::lexical_cast<Real>(to_species.get_attribute("D")));
    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.);
    for (std::vector<ReactionRule>::const_iterator itr(rules.begin());
            itr != rules.end(); ++itr)
    {
        const Real k((*itr).k());
        const Real P(k * factor/ (6 * sqrt(2.)
                    * (Da + Db) * world_->normalized_voxel_radius()));
        accp += P;
        if (accp >= rnd)
        {
            apply_reaction_(*itr, from_coord, to_coord);
        }
    }
}

void LatticeSimulator::apply_reaction_(const ReactionRule& reaction_rule,
        const Coord& from_coord, const Coord& to_coord)
{
    const ReactionRule::product_container_type& products(
            reaction_rule.products());
    if (products.size() == 0)
    {
        // Not yet test
        world_->remove_molecule(from_coord);
        world_->remove_molecule(to_coord);
    }
    else if (products.size() == 1)
    {
        Species product_species(*(products.begin()));
        world_->remove_molecule(from_coord);
        world_->update_molecule(to_coord, product_species);
    }
    else if (products.size() == 2)
    {
        // Not yet test
        const Integer rnd(world_->rng()->uniform_int(0,11));
        std::pair<Coord, bool> retval(world_->move_to_neighbor(to_coord, rnd));
        if (retval.second)
        {
            const Species product_species0(*(products.begin())),
                          product_species1(*(++(products.begin())));
            world_->remove_molecule(from_coord);
            world_->update_molecule(to_coord, product_species0);
            world_->update_molecule(retval.first, product_species1);
        }
    }
}

void LatticeSimulator::apply_reaction_(const ReactionRule& reaction_rule,
        const Coord& coord)
{
        if (reaction_rule.products().size() == 0)
        {
            world_->remove_molecule(coord);
        }
        else if (reaction_rule.products().size() == 1)
        {
            const Species product(*(reaction_rule.products().begin()));
            world_->update_molecule(coord, product);
        }
        else if (reaction_rule.products().size() == 2)
        {
            const Species product_species0(*(reaction_rule.products().begin())),
                          product_species1(*(++(reaction_rule.products().begin())));
            const Integer rnd(world_->rng()->uniform_int(0,11));
            std::pair<Coord, bool> retval(world_->move_to_neighbor(coord, rnd));
            if (retval.second)
            {
                world_->add_molecule(product_species0, coord);
                world_->update_molecule(retval.first, product_species1);
            }
        }

}

void LatticeSimulator::step()
{
    if (!is_initialized_)
    {
        initialize();
    }
    step_();
}

bool LatticeSimulator::step(const Real& upto)
{
    if (!is_initialized_)
    {
        initialize();
    }

    if (upto < t())
    {
        return false;
    }

    if (upto >= scheduler_.top().second->time())
    {
        step_();
        return true;
    }

    world_->set_t(upto);

    return false;
}

void LatticeSimulator::step_()
{
    EventScheduler::value_type top(scheduler_.pop());
    const Real time(top.second->time());
    top.second->fire(); // top.second->time_ is updated in fire()
    world_->set_t(time);
    EventScheduler::events_range events(scheduler_.events());
    for (EventScheduler::events_range::iterator itr(events.begin());
            itr != events.end(); ++itr)
    {
        (*itr).second->interrupt(time);
        scheduler_.update(*itr);
    }
    scheduler_.add(top.second);
}

} // lattice

} // ecell4
