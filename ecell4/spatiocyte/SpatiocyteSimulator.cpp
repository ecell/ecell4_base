#include "SpatiocyteSimulator.hpp"

#include <algorithm>
#include <iterator>
#include <ecell4/core/StructureType.hpp>

namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteSimulator::initialize()
{
    scheduler_.clear();
    update_alpha_map();
    const std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
        itr != species.end(); ++itr)
    {
        register_events(*itr);
    }


    const std::vector<ReactionRule>& rules(model_->reaction_rules());
    for (std::vector<ReactionRule>::const_iterator i(rules.begin());
        i != rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        if (rr.reactants().size() != 0)
        {
            continue;
        }
        const boost::shared_ptr<SpatiocyteEvent>
            zeroth_order_reaction_event(
                create_zeroth_order_reaction_event(rr, world_->t()));
        scheduler_.add(zeroth_order_reaction_event);
    }

    dt_ = scheduler_.next_time() - t();
}

void SpatiocyteSimulator::update_alpha_map()
{
    boost::shared_ptr<Model> model_(model());
    if (!model_ || !model_->is_static())
        return;

    const Model::reaction_rule_container_type reaction_rules(model_->reaction_rules());
    for (Model::reaction_rule_container_type::const_iterator itr(reaction_rules.begin());
            itr != reaction_rules.end(); ++itr)
    {
        const ReactionRule::reactant_container_type& reactants((*itr).reactants());
        if (reactants.size() != 2)
            continue;

        const Real alpha(calculate_alpha(*itr));
        for (int i(0); i < 2; ++i) {
            const Species& sp(reactants.at(i));
            alpha_map_type::iterator map_itr(alpha_map_.find(sp));
            if (map_itr == alpha_map_.end())
                alpha_map_.insert(alpha_map_type::value_type(sp, alpha));
            else if ((*map_itr).second > alpha)
                (*map_itr).second = alpha;
        }
    }
}

void SpatiocyteSimulator::register_events(const Species& sp)
{
    if (world_->has_molecule_pool(sp))
    {
        //TODO: Call steps only if sp is assigned not to StructureType.
        const boost::shared_ptr<SpatiocyteEvent> step_event(
                create_step_event(sp, world_->t()));
        scheduler_.add(step_event);
    }

    std::vector<ReactionRule> reaction_rules(model_->query_reaction_rules(sp));
    for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
        i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        const boost::shared_ptr<SpatiocyteEvent>
            first_order_reaction_event(
                create_first_order_reaction_event(rr, world_->t()));
        scheduler_.add(first_order_reaction_event);
    }
}

boost::shared_ptr<SpatiocyteEvent> SpatiocyteSimulator::create_step_event(
        const Species& species, const Real& t)
{
    double alpha(alpha_);
    alpha_map_type::const_iterator itr(alpha_map_.find(species));
    if (itr != alpha_map_.end() && (*itr).second < alpha)
        alpha = (*itr).second;

    boost::shared_ptr<SpatiocyteEvent> event(
        new StepEvent(this, model_, world_, species, t, alpha));
    return event;
}

boost::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_zeroth_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<SpatiocyteEvent> event(
            new ZerothOrderReactionEvent(world_, reaction_rule, t));
    return event;
}

boost::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_first_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<SpatiocyteEvent> event(new FirstOrderReactionEvent(
                world_, reaction_rule, t));
    return event;
}

void SpatiocyteSimulator::finalize()
{
    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator itr(events.begin());
            itr != events.end(); ++itr)
    {
        const Real queued_time((*itr).second->time() - (*itr).second->dt());
        StepEvent* step_event(dynamic_cast<StepEvent*>((*itr).second.get()));
        if (step_event != NULL && queued_time < t())
        {
            const Real alpha((t() - queued_time) / (*itr).second->dt());
            step_event->walk(alpha);
        }
    }
    initialize();
}

Real SpatiocyteSimulator::calculate_dimensional_factor(
    const VoxelPool* mt0, const VoxelPool* mt1, boost::shared_ptr<SpatiocyteWorld> world)
{
    const Species&
        speciesA(mt0->species()),
        speciesB(mt1->species());
    const Real
        D_A(mt0->D()),
        D_B(mt1->D());
    const Shape::dimension_kind
        dimensionA(mt0->get_dimension()),
        dimensionB(mt1->get_dimension());
    const Real Dtot(D_A + D_B);
    const Real gamma(pow(2 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0) + sqrt(22.0), 2) /
        (72 * (6 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0))));
    Real factor(0);
    if (dimensionA == Shape::THREE && dimensionB == Shape::THREE)
    {
        // if (speciesA != speciesB)
        //     factor = 1. / (6 * sqrt(2.0) * Dtot * world->voxel_radius());
        // else
        //     factor = 1. / (6 * sqrt(2.0) * D_A * world->voxel_radius());
        factor = 1. / (6 * sqrt(2.0) * Dtot * world->voxel_radius());
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::TWO)
    {
        // if (speciesA != speciesB)
        //     factor = gamma / Dtot;
        // else
        //     factor = gamma / D_A;
        factor = gamma / Dtot;
    }
    else if (dimensionA == Shape::THREE && dimensionB == Shape::TWO)
    {
        factor = sqrt(2.0) / (3 * D_A * world->voxel_radius());
        if (mt1->is_structure()) // B is Surface
        {
            factor *= world->unit_area();
        }
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::THREE)
    {
        factor = sqrt(2.0) / (3 * D_B * world->voxel_radius());
        if (mt0->is_structure()) // A is Surface
        {
            factor *= world->unit_area();
        }
    }
    else
        throw NotSupported("The dimension of a structure must be two or three.");
    return factor;
}

Real SpatiocyteSimulator::calculate_alpha(const ReactionRule& rule) const
{
    const ReactionRule::reactant_container_type& reactants(rule.reactants());
    if (reactants.size() != 2)
        return 1.0;

    const Species species[2] = {reactants.at(0), reactants.at(1)};
    const MoleculeInfo info[2] = {
        world_->get_molecule_info(species[0]),
        world_->get_molecule_info(species[1])
    };
    VoxelPool* mt[2];
    bool is_created[2] = {false, false};
    for (int i(0); i < 2; ++i) {
        try
        {
            mt[i] = world_->find_voxel_pool(species[i]);
        }
        catch(NotFound e)
        {
            VoxelPool *location(&(VacantType::getInstance()));
            if (info[i].loc != "") {
                try
                {
                    location = world_->find_voxel_pool(Species(info[i].loc));
                }
                catch(NotFound e)
                {
                    ;
                }
            }
            mt[i] = new MolecularType(species[i], location, info[i].radius, info[i].D);
            is_created[i] = true;
        }
    }
    const Real factor(calculate_dimensional_factor(mt[0], mt[1], world_));
    for (int i(0); i < 2; ++i)
        if (is_created[i])
            delete mt[i];
    const Real alpha(1.0 / (factor * rule.k()));
    return alpha < 1.0 ? alpha : 1.0;
}

/*
 * the Reaction between two molecules
 */
void SpatiocyteSimulator::step()
{
    step_();
    dt_ = scheduler_.next_time() - t();
}

bool SpatiocyteSimulator::step(const Real& upto)
{
    if (upto < t())
    {
        return false;
    }

    if (scheduler_.size() > 0 && upto >= scheduler_.top().second->time())
    {
        step_();
        dt_ = scheduler_.next_time() - t();
        return true;
    }

    world_->set_t(upto); //XXX: TODO
    last_reactions_.clear();
    dt_ = scheduler_.next_time() - t();
    return false;
}

void SpatiocyteSimulator::step_()
{
    last_reactions_.clear();

    scheduler_type::value_type top(scheduler_.pop());
    const Real time(top.second->time());
    world_->set_t(time);
    top.second->fire(); // top.second->time_ is updated in fire()
    if (!check_reaction())
        last_reactions_ = top.second->last_reactions();

    std::vector<Species> new_species;
    for (std::vector<reaction_type>::const_iterator itr(last_reactions_.begin());
            itr != last_reactions_.end(); ++itr)
        for (ReactionInfo::container_type::const_iterator
                product((*itr).second.products().begin());
                product != (*itr).second.products().end(); ++product)
        {
            const Species& species((*product).second.species());
            if (!world_->has_species(species))
                new_species.push_back(species);
        }

    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator itr(events.begin());
        itr != events.end(); ++itr)
    {
        (*itr).second->interrupt(time);
        scheduler_.update(*itr);
    }
    scheduler_.add(top.second);

    // update_alpha_map(); // may be performance cost
    for (std::vector<Species>::const_iterator itr(new_species.begin());
        itr != new_species.end(); ++itr)
    {
        register_events(*itr);
    }

    num_steps_++;
}

} // spatiocyte

} // ecell4
