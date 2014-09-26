#include "LatticeSimulator.hpp"

//#include <set>

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
        register_events(*itr);
    }

    // Model::reaction_rule_container_type rules(model_->reaction_rules());
    // for (Model::reaction_rule_container_type::iterator itr(rules.begin());
    //         itr != rules.end(); ++itr)
    // std::vector<ReactionRule> rules(model_->list_reaction_rules());
    // for (std::vector<ReactionRule>::iterator itr(rules.begin());
    //     itr != rules.end(); ++itr)
    // {
    //     if ((*itr).reactants().size() == 1)
    //     {
    //         const boost::shared_ptr<EventScheduler::Event> event(
    //                 create_first_order_reaction_event(*itr));
    //         scheduler_.add(event);
    //     }
    // }

    dt_ = scheduler_.next_time() - t();
}

void LatticeSimulator::finalize()
{
    EventScheduler::events_range events(scheduler_.events());
    for (EventScheduler::events_range::iterator itr(events.begin());
            itr != events.end(); ++itr)
    {
        const Real queued_time((*itr).second->time() - (*itr).second->dt());
        StepEvent* step_event(dynamic_cast<StepEvent*>((*itr).second.get()));
        if (step_event != NULL && queued_time < t())
        {
            const Real alpha((t() - queued_time) / (*itr).second->dt());
            walk(step_event->species(), alpha);
        }
    }
    initialize();
}

boost::shared_ptr<EventScheduler::Event> LatticeSimulator::create_step_event(
        const Species& species, const Real& t)
{
    boost::shared_ptr<EventScheduler::Event> event(new StepEvent(this, species, t));
    return event;
}

boost::shared_ptr<EventScheduler::Event>
LatticeSimulator::create_first_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<EventScheduler::Event> event(new FirstOrderReactionEvent(
                this, reaction_rule, t));
    return event;
}

std::pair<bool, Reaction<Voxel> > LatticeSimulator::attempt_reaction_(
    const LatticeWorld::particle_info info, LatticeWorld::coordinate_type to_coord)
{
    const Species
        from_species(world_->get_molecular_type_private(info.first)->species());
    const MolecularTypeBase* to_mt(world_->get_molecular_type_private(to_coord));
    const Species to_species(to_mt->species());
    const std::vector<ReactionRule> rules(model_->query_reaction_rules(
                from_species, to_species));

    const Real factor((from_species == to_species) ? 2 : 1);
    const LatticeWorld::molecule_info_type
        from_minfo(world_->get_molecule_info(from_species)),
        to_minfo(world_->get_molecule_info(to_species));

    const Real Dtot(from_minfo.D + to_minfo.D);
    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.);
    for (std::vector<ReactionRule>::const_iterator itr(rules.begin());
            itr != rules.end(); ++itr)
    {
        const Real k((*itr).k());
        const Real P(k * factor/ (6 * sqrt(2.)
                    * Dtot * world_->voxel_radius()));
        accp += P;
        if (accp > 1)
        {
            std::cerr << "accp : " << accp << std::endl;
        }
        if (accp >= rnd)
        {
            LatticeWorld::particle_info to_info(*(to_mt->find(to_coord)));
            return apply_reaction_(*itr, info, to_info);
        }
    }
    return std::pair<bool, Reaction<Voxel> >(false, Reaction<Voxel>());
}

/*
 * the Reaction between two molecules
 */
std::pair<bool, Reaction<Voxel> > LatticeSimulator::apply_reaction_(
    const ReactionRule& reaction_rule,
    const LatticeWorld::particle_info from_info,
    const LatticeWorld::particle_info to_info)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());
    Reaction<Voxel> reaction;
    reaction.rule = reaction_rule;
    MolecularTypeBase* from_mtype(world_->get_molecular_type_private(from_info.first));
    MolecularTypeBase* to_mtype(world_->get_molecular_type_private(to_info.first));

    const std::string from_loc((from_mtype->location()->is_vacant())
        ? "" : from_mtype->location()->species().serial());
    const std::string to_loc((to_mtype->location()->is_vacant())
        ? "" : to_mtype->location()->species().serial());

    reaction.reactants.push_back(
        Reaction<Voxel>::particle_type(
            from_info.second,
            Voxel(from_mtype->species(), world_->private2coord(from_info.first),
                from_mtype->radius(), from_mtype->D(), from_loc)));
    reaction.reactants.push_back(
        Reaction<Voxel>::particle_type(
            to_info.second,
            Voxel(to_mtype->species(), world_->private2coord(to_info.first),
                to_mtype->radius(), to_mtype->D(), to_loc)));

    if (products.size() == 0)
    {
        // Not test yet
        world_->remove_voxel_private(from_info.first);
        world_->remove_voxel_private(to_info.first);
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (products.size() == 1)
    {
        const LatticeWorld::private_coordinate_type coord(to_info.first);
        Species product_species(*(products.begin()));
        world_->remove_voxel_private(from_info.first);
        world_->remove_voxel_private(to_info.first);

        if (!world_->has_species(product_species))
            new_species_.push_back(product_species);

        const LatticeWorld::molecule_info_type
            minfo(world_->get_molecule_info(product_species));
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel_private(
                Voxel(product_species, coord, minfo.radius, minfo.D, minfo.loc)));
        if (new_mol.second)
        {
            return std::pair<bool, Reaction<Voxel> >(false, reaction);
        }

        reaction.products.push_back(
            Reaction<Voxel>::particle_type(
                new_mol.first.first,
                Voxel(product_species, world_->private2coord(coord),
                    minfo.radius, minfo.D, minfo.loc)));
        reactions_.push_back(reaction_rule);
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (products.size() == 2)
    {
        // Not test yet
        const LatticeWorld::private_coordinate_type from_coord(from_info.first);
        const LatticeWorld::private_coordinate_type to_coord(to_info.first);
        world_->remove_voxel_private(from_coord);
        world_->remove_voxel_private(to_coord);

        const Species product_species0(*(products.begin())),
                      product_species1(*(++(products.begin())));

        if (!world_->has_species(product_species0))
            new_species_.push_back(product_species0);
        const LatticeWorld::molecule_info_type
            minfo0(world_->get_molecule_info(product_species0));
        std::pair<std::pair<ParticleID, Voxel>, bool> new0(
            world_->new_voxel_private(
                Voxel(product_species0, from_coord, minfo0.radius, minfo0.D, minfo0.loc)));
        reaction.products.push_back(
            Reaction<Voxel>::particle_type(
                new0.first.first,
                Voxel(product_species0, world_->private2coord(from_coord),
                    minfo0.radius, minfo0.D, minfo0.loc)));

        if (!world_->has_species(product_species1))
            new_species_.push_back(product_species1);
        const LatticeWorld::molecule_info_type
            minfo1(world_->get_molecule_info(product_species1));
        std::pair<std::pair<ParticleID, Voxel>, bool> new1(
            world_->new_voxel_private(
                Voxel(product_species1, to_coord, minfo1.radius, minfo1.D, minfo1.loc)));
        reaction.products.push_back(
            Reaction<Voxel>::particle_type(
                new1.first.first,
                Voxel(product_species1, world_->private2coord(to_coord),
                    minfo1.radius, minfo1.D, minfo1.loc)));

        reactions_.push_back(reaction_rule);
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    return std::pair<bool, Reaction<Voxel> >(false, reaction);
}

/*
 * the First Order Reaction
 */
std::pair<bool, Reaction<Voxel> > LatticeSimulator::apply_reaction_(
        const ReactionRule& reaction_rule, const LatticeWorld::particle_info info)
{
    Reaction<Voxel> reaction;
    reaction.rule = reaction_rule;
    MolecularTypeBase* mtype(world_->get_molecular_type_private(info.first));
    const std::string src_loc((mtype->location()->is_vacant())
        ? "" : mtype->location()->species().serial());
    reaction.reactants.push_back(
        Reaction<Voxel>::particle_type(
            info.second,
            Voxel(mtype->species(), world_->private2coord(info.first),
                mtype->radius(), mtype->D(), src_loc)));

    if (reaction_rule.products().size() == 0)
    {
        world_->remove_voxel_private(info.first);
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (reaction_rule.products().size() == 1)
    {
        const Species product(*(reaction_rule.products().begin()));
        if (!world_->has_species(product))
                new_species_.push_back(product);
        const LatticeWorld::molecule_info_type
            minfo(world_->get_molecule_info(product));
        const LatticeWorld::private_coordinate_type coord(info.first);
        world_->remove_voxel_private(coord);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_val(
            world_->new_voxel_private(
                Voxel(product, coord, minfo.radius, minfo.D, minfo.loc)));

        reactions_.push_back(reaction_rule);

        reaction.products.push_back(
            std::pair<ParticleID, Voxel>(
                new_val.first.first,
                Voxel(product, world_->private2coord(coord),
                    minfo.radius, minfo.D, minfo.loc)));
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (reaction_rule.products().size() == 2)
    {
        const Species product_species0(*(reaction_rule.products().begin())),
                      product_species1(*(++(reaction_rule.products().begin())));
        const LatticeWorld::private_coordinate_type coord(info.first);
        std::pair<LatticeWorld::private_coordinate_type, bool> neighbor(
                world_->check_neighbor_private(coord));
        if (neighbor.second)
        {
            if (!world_->has_species(product_species0))
                new_species_.push_back(product_species0);
            const LatticeWorld::molecule_info_type
                minfo0(world_->get_molecule_info(product_species0));
            world_->remove_voxel_private(coord);
            std::pair<std::pair<ParticleID, Voxel>, bool> new_val0(
                world_->new_voxel_private(
                    Voxel(product_species0, coord, minfo0.radius, minfo0.D, minfo0.loc)));

            if (!world_->has_species(product_species1))
                new_species_.push_back(product_species1);
            const LatticeWorld::molecule_info_type
                minfo1(world_->get_molecule_info(product_species1));
            std::pair<std::pair<ParticleID, Voxel>, bool> new_val1(
                world_->new_voxel_private(
                    Voxel(product_species1, neighbor.first, minfo1.radius, minfo1.D,
                        minfo1.loc)));

            reactions_.push_back(reaction_rule);

            reaction.products.push_back(
                std::pair<ParticleID, Voxel>(
                    new_val0.first.first,
                    Voxel(product_species0, world_->private2coord(coord),
                        minfo0.radius, minfo0.D, minfo0.loc)));
            reaction.products.push_back(
                std::pair<ParticleID, Voxel>(
                    new_val1.first.first,
                    Voxel(product_species1, world_->private2coord(neighbor.first),
                        minfo1.radius, minfo1.D, minfo1.loc)));
            return std::pair<bool, Reaction<Voxel> >(true, reaction);
        }
    }
    return std::pair<bool, Reaction<Voxel> >(false, reaction);
}

// void LatticeSimulator::register_step_event(const Species& species)
// {
//     const boost::shared_ptr<EventScheduler::Event> event(
//             create_step_event(species, world_->t()));
//     scheduler_.add(event);
// }

void LatticeSimulator::register_events(const Species& sp)
{
    const boost::shared_ptr<EventScheduler::Event> step_event(
            create_step_event(sp, world_->t()));
    scheduler_.add(step_event);

    std::vector<ReactionRule> reaction_rules(model_->query_reaction_rules(sp));
    for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
        i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        const boost::shared_ptr<EventScheduler::Event>
            first_order_reaction_event(
                create_first_order_reaction_event(rr, world_->t()));
        scheduler_.add(first_order_reaction_event);
    }
}

void LatticeSimulator::step()
{
    step_();
    dt_ = scheduler_.next_time() - t();
}

bool LatticeSimulator::step(const Real& upto)
{
    if (upto < t())
    {
        return false;
    }

    if (upto >= scheduler_.top().second->time())
    {
        step_();
        dt_ = scheduler_.next_time() - t();
        return true;
    }

    world_->set_t(upto); //XXX: TODO
    reactions_.clear();
    new_species_.clear();
    dt_ = scheduler_.next_time() - t();
    return false;
}

void LatticeSimulator::step_()
{
    reactions_.clear();
    new_species_.clear();

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

    for (std::vector<Species>::const_iterator itr(new_species_.begin());
            itr != new_species_.end(); ++itr)
    {
        register_events(*itr);
    }

    num_steps_++;
}

// void LatticeSimulator::run(const Real& duration)
// {
//     initialize();
//     const Real upto(t() + duration);
//     while (step(upto))
//         ;
//     finalize();
// }

void LatticeSimulator::walk(const Species& species)
{
    walk(species, 1.0);
}

void LatticeSimulator::walk(const Species& species, const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
        return; // INVALID ALPHA VALUE

    boost::shared_ptr<RandomNumberGenerator> rng(world_->rng());

    MolecularTypeBase* mtype(world_->find_molecular_type(species));
    mtype->shuffle(*rng);
    std::vector<ParticleID> pids;
    Integer i(0);
    Integer max(rng->binomial(alpha, mtype->size()));
    pids.reserve(max);
    while(i < max)
    {
        const std::pair<std::pair<LatticeWorld::particle_info,
            LatticeWorld::private_coordinate_type>, bool> neighbor(
                    world_->move_to_neighbor(mtype, i));
        const LatticeWorld::particle_info info(neighbor.first.first);
        const LatticeWorld::private_coordinate_type coord(neighbor.first.second);
        pids.push_back(info.second);

        if (!neighbor.second)
        {
            const std::pair<bool, Reaction<Voxel> > retval(attempt_reaction_(info, coord));
            if (retval.first)
            {
                const Reaction<Voxel> reaction(retval.second);
                for (std::vector<Reaction<Voxel>::particle_type>::const_iterator itr(reaction.reactants.begin());
                        itr != reaction.reactants.end(); ++itr)
                {
                    for (std::vector<ParticleID>::const_iterator pid_itr(pids.begin());
                            pid_itr != pids.end(); ++pid_itr)
                    {
                        if (*pid_itr == (*itr).first)
                        {
                            --i;
                            break;
                        }
                    }
                }
            }
        }
        ++i;
        if (max > mtype->size())
            max = mtype->size();
    }
}

} // lattice

} // ecell4
