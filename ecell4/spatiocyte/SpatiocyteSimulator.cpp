#include "SpatiocyteSimulator.hpp"
#include "utils.hpp"

#include <algorithm>
#include <ecell4/core/StructureType.hpp>
#include <iterator>

namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteSimulator::initialize()
{
    last_reactions_.clear();
    species_list_.clear(); // XXX:FIXME: Messy patch

    scheduler_.clear();
    update_alpha_map();
    for (const auto &species : world_->list_species())
    {
        register_events(species);
    }

    for (const auto &rule : model_->reaction_rules())
    {
        if (rule.reactants().size() != 0)
        {
            continue;
        }
        const std::shared_ptr<SpatiocyteEvent> zeroth_order_reaction_event(
            create_zeroth_order_reaction_event(rule, world_->t()));
        scheduler_.add(zeroth_order_reaction_event);
    }

    dt_ = scheduler_.next_time() - t();
}

void SpatiocyteSimulator::update_alpha_map()
{
    std::shared_ptr<Model> model_(model());
    if (!model_ || !model_->is_static())
        return;

    for (const auto &rule : model_->reaction_rules())
    {
        const ReactionRule::reactant_container_type &reactants(
            rule.reactants());
        if (reactants.size() != 2)
            continue;

        const Real alpha(calculate_alpha(rule, world_));
        for (int i(0); i < 2; ++i)
        {
            const Species &sp(reactants.at(i));
            alpha_map_type::iterator map_itr(alpha_map_.find(sp));
            if (map_itr == alpha_map_.end())
                alpha_map_.insert(alpha_map_type::value_type(sp, alpha));
            else if ((*map_itr).second > alpha)
                (*map_itr).second = alpha;
        }
    }
}

void SpatiocyteSimulator::register_events(const Species &sp)
{
    species_list_.push_back(sp); // XXX:FIXME: Messy patch

    if (world_->has_molecule_pool(sp))
    {
        // TODO: Call steps only if sp is assigned not to StructureType.
        alpha_map_type::const_iterator itr(alpha_map_.find(sp));
        const Real alpha(itr != alpha_map_.end() ? itr->second : 1.0);
        const std::shared_ptr<SpatiocyteEvent> step_event(
            create_step_event(sp, world_->t(), alpha));
        scheduler_.add(step_event);
    }

    for (const auto &rule : model_->query_reaction_rules(sp))
    {
        const std::shared_ptr<SpatiocyteEvent> first_order_reaction_event(
            create_first_order_reaction_event(rule, world_->t()));
        scheduler_.add(first_order_reaction_event);
    }
}

std::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_step_event(const Species &species, const Real &t,
                                       const Real &alpha)
{
    std::shared_ptr<MoleculePool> mpool(world_->find_molecule_pool(species));
    const Shape::dimension_kind dimension(world_->get_dimension(species));

    if (dimension == Shape::THREE)
    {
        return std::shared_ptr<SpatiocyteEvent>(
            new StepEvent<3>(model_, world_, species, t, alpha));
    }
    else if (dimension == Shape::TWO)
    {
        return std::shared_ptr<SpatiocyteEvent>(
            new StepEvent<2>(model_, world_, species, t, alpha));
    }
    else
    {
        throw NotSupported(
            "The dimension of a structure must be two or three.");
    }
}

std::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_zeroth_order_reaction_event(
    const ReactionRule &reaction_rule, const Real &t)
{
    std::shared_ptr<SpatiocyteEvent> event(
        new ZerothOrderReactionEvent(world_, reaction_rule, t));
    return event;
}

std::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_first_order_reaction_event(
    const ReactionRule &reaction_rule, const Real &t)
{
    std::shared_ptr<SpatiocyteEvent> event(
        new FirstOrderReactionEvent(world_, reaction_rule, t));
    return event;
}

void SpatiocyteSimulator::finalize()
{
    for (const auto &item : scheduler_.events())
    {
        const auto &event(item.second);
        const Real queued_time(event->time() - event->dt());
        auto *step3_event(dynamic_cast<StepEvent<3> *>(event.get()));
        if (step3_event != NULL && queued_time < t())
        {
            const Real factor((t() - queued_time) / event->dt());
            // assert(factor <= 1);
            step3_event->walk(step3_event->alpha() * factor);
        }
        auto *step2_event(dynamic_cast<StepEvent<2> *>(event.get()));
        if (step2_event != NULL && queued_time < t())
        {
            const Real factor((t() - queued_time) / event->dt());
            // assert(factor <= 1);
            step2_event->walk(step2_event->alpha() * factor);
        }
    }

    initialize();
}

void SpatiocyteSimulator::step()
{
    step_();
    dt_ = scheduler_.next_time() - t();
}

bool SpatiocyteSimulator::step(const Real &upto)
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

    world_->set_t(upto);
    last_reactions_.clear();
    dt_ = scheduler_.next_time() - t();
    finalize();
    return false;
}

void SpatiocyteSimulator::step_()
{
    scheduler_type::value_type top(scheduler_.pop());
    const Real time(top.second->time());
    world_->set_t(time);
    top.second->fire(); // top.second->time_ is updated in fire()
    set_last_event_(
        std::const_pointer_cast<const SpatiocyteEvent>(top.second));

    last_reactions_ = last_event_->reactions();

    std::vector<Species> new_species;
    for (const auto &reaction : last_reactions())
    {
        for (const auto &product : reaction.second.products())
        {
            const Species &species(product.species);
            // if (!world_->has_species(species))
            if (std::find(species_list_.begin(), species_list_.end(),
                          species) ==
                species_list_.end()) // XXX:FIXME: Messy patch
                new_species.push_back(species);
        }
    }

    for (const auto &event : scheduler_.events())
    {
        event.second->interrupt(time);
        scheduler_.update(event);
    }
    scheduler_.add(top.second);

    // update_alpha_map(); // may be performance cost
    for (const auto &species : new_species)
    {
        register_events(species);
    }

    num_steps_++;
}

} // namespace spatiocyte

} // namespace ecell4
