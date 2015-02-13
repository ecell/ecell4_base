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

std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::attempt_reaction_(
    const LatticeWorld::particle_info_type info, LatticeWorld::coordinate_type to_coord)
{
    const MolecularTypeBase* from_mt(
        world_->get_molecular_type_private(info.first));
    const MolecularTypeBase* to_mt(
        world_->get_molecular_type_private(to_coord));

    if (to_mt->is_vacant())
    {
        return std::pair<bool, reaction_type>(false, reaction_type());
    }

    const Species
        from_species(from_mt->species()),
        to_species(to_mt->species());
    const LatticeWorld::molecule_info_type
        from_minfo(world_->get_molecule_info(from_mt)),
        to_minfo(world_->get_molecule_info(to_mt));
    // const LatticeWorld::molecule_info_type
    //     from_minfo(world_->get_molecule_info(from_species)),
    //     to_minfo(world_->get_molecule_info(to_species));

    const std::vector<ReactionRule> rules(
        model_->query_reaction_rules(from_species, to_species));

    const Real D0(from_minfo.D);
    const Real D1(to_minfo.D);
    const Shape::dimension_kind dimensionA(world_->get_dimension_kind(from_minfo.loc));
    const Shape::dimension_kind dimensionB(world_->get_dimension_kind(to_minfo.loc));

    const Real Dtot(D0 + D1);
    const Real rnd(world_->rng()->uniform(0,1));
    const Real gamma(pow(2*sqrt(2) + 4*sqrt(3) + 3*sqrt(6) + sqrt(22), 2) /
        (72*(6*sqrt(2) + 4*sqrt(3) + 3*sqrt(6))));
    Real factor(0);
    if (dimensionA == Shape::THREE && dimensionB == Shape::THREE)
    {
        if (from_species != to_species)
            factor = 1. / (6 * sqrt(2.) * Dtot * world_->voxel_radius());
        else
            factor = 1. / (6 * sqrt(2.) * D0 * world_->voxel_radius());
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::TWO)
    {
        if (from_species != to_species)
            factor = gamma / Dtot;
        else
            factor = gamma / D0;
    }
    else if (dimensionA == Shape::THREE && dimensionB == Shape::TWO)
    {
        factor = sqrt(2) / (3 * D0 * world_->voxel_radius());
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::THREE)
    {
        factor = sqrt(2) / (3 * D1 * world_->voxel_radius()); // 不要?
    }
    else
        throw NotSupported("The dimension of a shape must be two or three.");

    Real accp(0.0);
    for (std::vector<ReactionRule>::const_iterator itr(rules.begin());
        itr != rules.end(); ++itr)
    {
        const Real k((*itr).k());
        const Real P(k * factor);
        accp += P;
        if (accp > 1)
        {
            std::cerr << "accp : " << accp << std::endl;
        }
        if (accp >= rnd)
        {
            return apply_second_order_reaction_(*itr,
                world_->make_pid_voxel_pair(from_mt, info),
                world_->make_pid_voxel_pair(to_mt, to_coord));
        }
    }
    return std::pair<bool, reaction_type>(false, reaction_type());
}

/*
 * the Reaction between two molecules
 */
std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::apply_second_order_reaction_(
    const ReactionRule& reaction_rule,
    const LatticeSimulator::reaction_type::particle_type& p0,
    const LatticeSimulator::reaction_type::particle_type& p1)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;
    reaction.reactants.push_back(p0);
    reaction.reactants.push_back(p1);

    const LatticeWorld::private_coordinate_type from_coord(
        world_->coord2private(p0.second.coordinate()));
    const LatticeWorld::private_coordinate_type to_coord(
        world_->coord2private(p1.second.coordinate()));

    world_->remove_voxel_private(from_coord);
    world_->remove_voxel_private(to_coord);

    const Species& product_species0(*(products.begin()));
    const Species& product_species1(*(++(products.begin())));
    switch(products.size())
    {
        case 1:
            apply_ab2c(to_coord, product_species0, reaction);
            break;
        case 2:
            apply_ab2cd(from_coord, to_coord, product_species0, product_species1, reaction);
            break;
        default:
            return std::pair<bool, reaction_type>(false, reaction);
    }

    reactions_.push_back(reaction_rule);
    // for (ReactionRule::product_container_type::const_iterator
    //     i(products.begin()); i != products.end(); ++i)
    // {
    //     register_product_species(*i);
    // }

    return std::pair<bool, reaction_type>(true, reaction);
}

// std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::apply_second_order_reaction_(
//     const ReactionRule& reaction_rule,
//     const LatticeWorld::particle_info_type from_info,
//     const LatticeWorld::particle_info_type to_info)
// {
//     const ReactionRule::product_container_type&
//         products(reaction_rule.products());
//     reaction_type reaction;
//     reaction.rule = reaction_rule;
// 
//     MolecularTypeBase* from_mtype(world_->get_molecular_type_private(from_info.first));
//     MolecularTypeBase* to_mtype(world_->get_molecular_type_private(to_info.first));
// 
//     const std::string from_loc((from_mtype->location()->is_vacant())
//         ? "" : from_mtype->location()->species().serial());
//     const std::string to_loc((to_mtype->location()->is_vacant())
//         ? "" : to_mtype->location()->species().serial());
// 
//     reaction.reactants.push_back(
//         reaction_type::particle_type(
//             from_info.second,
//             Voxel(from_mtype->species(), world_->private2coord(from_info.first),
//                 from_mtype->radius(), from_mtype->D(), from_loc)));
//     reaction.reactants.push_back(
//         reaction_type::particle_type(
//             to_info.second,
//             Voxel(to_mtype->species(), world_->private2coord(to_info.first),
//                 to_mtype->radius(), to_mtype->D(), to_loc)));
// 
//     const LatticeWorld::private_coordinate_type from_coord(from_info.first);
//     const LatticeWorld::private_coordinate_type to_coord(to_info.first);
// 
//     world_->remove_voxel_private(from_coord);
//     world_->remove_voxel_private(to_coord);
// 
//     const Species& product_species0(*(products.begin()));
//     const Species& product_species1(*(++(products.begin())));
//     switch(products.size())
//     {
//         case 1:
//             apply_ab2c(to_coord, product_species0, reaction);
//             break;
//         case 2:
//             apply_ab2cd(from_coord, to_coord, product_species0, product_species1, reaction);
//             break;
//         default:
//             return std::pair<bool, reaction_type>(false, reaction);
//     }
// 
//     reactions_.push_back(reaction_rule);
//     // for (ReactionRule::product_container_type::const_iterator
//     //     i(products.begin()); i != products.end(); ++i)
//     // {
//     //     register_product_species(*i);
//     // }
// 
//     return std::pair<bool, reaction_type>(true, reaction);
// }

void LatticeSimulator::apply_ab2c(
    const LatticeWorld::private_coordinate_type coord,
    const Species& product_species,
    reaction_type& reaction)
{
    register_product_species(product_species);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
        world_->new_voxel_private(product_species, coord));
    if (!new_mol.second)
    {
        throw IllegalState("no place for the product.");
    }
    reaction.products.push_back(
        reaction_type::particle_type(
            new_mol.first.first,
            this->private_voxel2voxel(new_mol.first.second)));
}

// Not tested yet
void LatticeSimulator::apply_ab2cd(
    const LatticeWorld::private_coordinate_type from_coord,
    const LatticeWorld::private_coordinate_type to_coord,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    register_product_species(product_species0);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
        world_->new_voxel_private(product_species0, from_coord));
    if (!new_mol0.second)
    {
        throw IllegalState("no place for the first product.");
    }
    reaction.products.push_back(
        reaction_type::particle_type(
            new_mol0.first.first,
            this->private_voxel2voxel(new_mol0.first.second)));

    register_product_species(product_species1);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
        world_->new_voxel_private(product_species1, to_coord));
    if (!new_mol1.second)
    {
        throw IllegalState("no place for the second product.");
    }
    reaction.products.push_back(
        reaction_type::particle_type(
            new_mol1.first.first,
            this->private_voxel2voxel(new_mol1.first.second)));
}

/*
 * the First Order Reaction
 */
std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::apply_first_order_reaction_(
    const ReactionRule& reaction_rule,
    const LatticeSimulator::reaction_type::particle_type& p)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;
    reaction.reactants.push_back(p);

    const LatticeWorld::private_coordinate_type coord(
        world_->coord2private(p.second.coordinate()));
    const Species& product_species0(*(products.begin()));
    const Species& product_species1(*(++(products.begin())));
    switch(products.size()) {
        case 0:
            world_->remove_voxel_private(coord);
            break;
        case 1:
            apply_a2b(coord, product_species0, reaction);
            break;
        case 2:
            if (!apply_a2bc(coord, product_species0, product_species1, reaction)) {
                return std::pair<bool, reaction_type>(false, reaction);
            }
            break;
        default:
            return std::pair<bool, reaction_type>(false, reaction);
    }

    reactions_.push_back(reaction_rule);
    // for (ReactionRule::product_container_type::const_iterator
    //     i(products.begin()); i != products.end(); ++i)
    // {
    //     register_product_species(*i);
    // }

    return std::pair<bool, reaction_type>(true, reaction);
}

// std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::apply_first_order_reaction_(
//         const ReactionRule& reaction_rule, const LatticeWorld::particle_info_type info)
// {
//     const ReactionRule::product_container_type&
//         products(reaction_rule.products());
//     reaction_type reaction;
//     reaction.rule = reaction_rule;
//     MolecularTypeBase* mtype(world_->get_molecular_type_private(info.first));
//     const std::string src_loc((mtype->location()->is_vacant())
//         ? "" : mtype->location()->species().serial());
//     reaction.reactants.push_back(
//         reaction_type::particle_type(
//             info.second,
//             Voxel(mtype->species(), world_->private2coord(info.first),
//                 mtype->radius(), mtype->D(), src_loc)));
// 
//     const LatticeWorld::private_coordinate_type coord(info.first);
//     const Species& product_species0(*(products.begin()));
//     const Species& product_species1(*(++(products.begin())));
//     switch(products.size()) {
//         case 0:
//             world_->remove_voxel_private(info.first);
//             break;
//         case 1:
//             apply_a2b(coord, product_species0, reaction);
//             break;
//         case 2:
//             if (!apply_a2bc(coord, product_species0, product_species1, reaction)) {
//                 return std::pair<bool, reaction_type>(false, reaction);
//             }
//             break;
//         default:
//             return std::pair<bool, reaction_type>(false, reaction);
//     }
// 
//     reactions_.push_back(reaction_rule);
//     // for (ReactionRule::product_container_type::const_iterator
//     //     i(products.begin()); i != products.end(); ++i)
//     // {
//     //     register_product_species(*i);
//     // }
// 
//     return std::pair<bool, reaction_type>(true, reaction);
// }

void LatticeSimulator::apply_a2b(
    const LatticeWorld::private_coordinate_type coord,
    const Species& product_species,
    reaction_type& reaction)
{
    world_->remove_voxel_private(coord);

    register_product_species(product_species);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
        world_->new_voxel_private(product_species, coord));
    if (!new_mol.second)
    {
        throw IllegalState("no place for the product.");
    }
    reaction.products.push_back(
        reaction_type::particle_type(
            new_mol.first.first,
            this->private_voxel2voxel(new_mol.first.second)));
}

bool LatticeSimulator::apply_a2bc(
    const LatticeWorld::private_coordinate_type coord,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    std::pair<LatticeWorld::private_coordinate_type, bool> neighbor(
            world_->check_neighbor_private(coord));
    if (!neighbor.second)
    {
        return false;
    }

    world_->remove_voxel_private(coord);

    register_product_species(product_species0);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
        world_->new_voxel_private(product_species0, coord));
    if (!new_mol0.second)
    {
        throw IllegalState("no place for the first product.");
    }
    reaction.products.push_back(
        reaction_type::particle_type(
            new_mol0.first.first,
            this->private_voxel2voxel(new_mol0.first.second)));

    register_product_species(product_species1);
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
        world_->new_voxel_private(product_species1, neighbor.first));
    if (!new_mol1.second)
    {
        throw IllegalState("no place for the second product.");
    }
    reaction.products.push_back(
        reaction_type::particle_type(
            new_mol1.first.first,
            this->private_voxel2voxel(new_mol1.first.second)));
    return true;
}

void LatticeSimulator::register_product_species(const Species& product_species)
{
    if (!world_->has_species(product_species))
    {
        new_species_.push_back(product_species);
    }
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
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());

    MolecularTypeBase* mtype(world_->find_molecular_type(species));
    MolecularTypeBase* loc(mtype->location());
    //XXX: mtype->shuffle(*rng);

    Integer i(0), max(rng->binomial(alpha, mtype->size()));
    while (i < max)
    {
        const Integer rnd(rng->uniform_int(0, 11));
        const std::pair<LatticeWorld::private_coordinate_type, bool>
            neighbor(world_->move_to_neighbor(
                mtype, loc, (*mtype)[i], rnd));

        if (!neighbor.second)
        {
            const LatticeWorld::particle_info_type info((*mtype)[i]);
            const LatticeWorld::private_coordinate_type to_coord(neighbor.first);

            const std::pair<bool, reaction_type>
                retval(attempt_reaction_(info, to_coord));
            if (retval.first)
            {
                --i;
                --max;

                if (max > mtype->size())
                {
                    max = mtype->size(); //XXX: for a dimerization
                }
            }
        }

        ++i;
    }
}

} // lattice

} // ecell4
