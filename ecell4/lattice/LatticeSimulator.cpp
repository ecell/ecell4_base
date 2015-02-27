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
    const LatticeWorld::particle_info info, LatticeWorld::coordinate_type to_coord)
{
    const Species
        from_species(world_->get_molecular_type_private(info.first)->species());
    const MolecularTypeBase* to_mt(world_->get_molecular_type_private(to_coord));
    const Species to_species(to_mt->species());
    const std::vector<ReactionRule> rules(model_->query_reaction_rules(
                from_species, to_species));

    const LatticeWorld::molecule_info_type
        from_minfo(world_->get_molecule_info(from_species)),
        to_minfo(world_->get_molecule_info(to_species));

    const Real Dtot(from_minfo.D + to_minfo.D);
    const Real rnd(world_->rng()->uniform(0,1));
    const Shape::dimension_kind dimensionA(world_->get_dimension_kind(from_minfo.loc));
    const Shape::dimension_kind dimensionB(world_->get_dimension_kind(to_minfo.loc));
    const Real gamma(pow(2*sqrt(2) + 4*sqrt(3) + 3*sqrt(6) + sqrt(22), 2) /
        (72*(6*sqrt(2) + 4*sqrt(3) + 3*sqrt(6))));
    Real factor(0);
    if (dimensionA == Shape::THREE && dimensionB == Shape::THREE)
    {
        if (from_species != to_species)
            factor = 1. / (6 * sqrt(2.) * Dtot * world_->voxel_radius());
        else
            factor = 1. / (6 * sqrt(2.) * from_minfo.D * world_->voxel_radius());
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::TWO)
    {
        if (from_species != to_species)
            factor = gamma / Dtot;
        else
            factor = gamma / from_minfo.D;
    }
    else if (dimensionA == Shape::THREE && dimensionB == Shape::TWO)
    {
        factor = sqrt(2) / (3 * from_minfo.D * world_->voxel_radius());
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::THREE)
    {
        factor = sqrt(2) / (3 * to_minfo.D * world_->voxel_radius()); // 不要?
    }
    else
        throw NotSupported("The dimension of a shape must be two or three.");

    Real accp(0.);
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
            LatticeWorld::particle_info to_info(*(to_mt->find(to_coord)));
            return apply_second_order_reaction_(*itr, info, to_info);
        }
    }
    return std::pair<bool, reaction_type>(false, reaction_type());
}

/*
 * the Reaction between two molecules
 */
std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::apply_second_order_reaction_(
    const ReactionRule& reaction_rule,
    const LatticeWorld::particle_info from_info,
    const LatticeWorld::particle_info to_info)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;
    const MolecularTypeBase* from_mtype(
            world_->get_molecular_type_private(from_info.first));
    const MolecularTypeBase* to_mtype(
            world_->get_molecular_type_private(to_info.first));

    const LatticeWorld::private_coordinate_type from_coord(from_info.first);
    const LatticeWorld::private_coordinate_type to_coord(to_info.first);
    const Species& product_species0(*(products.begin()));
    const Species& product_species1(*(++(products.begin())));
    switch(products.size())
    {
        case 1:
            apply_ab2c(from_info, to_info, product_species0, reaction);
            break;
        case 2:
            apply_ab2cd(from_info, to_info,
                    product_species0, product_species1, reaction);
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

void LatticeSimulator::apply_ab2c(
    const LatticeWorld::particle_info from_info,
    const LatticeWorld::particle_info to_info,
    const Species& product_species,
    reaction_type& reaction)
{
    const std::string toloc(world_->get_molecule_info(product_species).loc);
    const std::string fserial(get_serial(from_info.first));
    const std::string tserial(get_serial(to_info.first));

    if (fserial == toloc)
    {
        register_reactant_species(from_info, reaction);
        register_reactant_species(to_info, reaction);
        register_product_species(product_species);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel_private(product_species, from_info.first));
        if (!new_mol.second)
        {
            throw IllegalState("no place for the product. 1");
        }
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol.first.first,
                this->private_voxel2voxel(new_mol.first.second)));
    }
    else if(tserial == toloc)
    {
        register_reactant_species(from_info, reaction);
        register_reactant_species(to_info, reaction);
        register_product_species(product_species);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel_private(product_species, to_info.first));
        if (!new_mol.second)
        {
            throw IllegalState("no place for the product. 2");
        }
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol.first.first,
                this->private_voxel2voxel(new_mol.first.second)));
    }
    else
    {
        const std::string state(
                "fserial: " + fserial +
                ", tserial: " + tserial +
                ", location: " + toloc);
        // throw IllegalState("no place for the product. 3");
        throw IllegalState(state);
    }
}

// Not tested yet
void LatticeSimulator::apply_ab2cd(
    const LatticeWorld::particle_info from_info,
    const LatticeWorld::particle_info to_info,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    const LatticeWorld::private_coordinate_type from_coord(from_info.first);
    const LatticeWorld::private_coordinate_type to_coord(to_info.first);
    const std::string aserial(get_serial(from_coord));
    const std::string bserial(get_serial(to_coord));
    const std::string cloc(world_->get_molecule_info(product_species0).loc);
    const std::string dloc(world_->get_molecule_info(product_species1).loc);

    if (aserial == cloc)
    {
        if (bserial == dloc)
        {
            register_reactant_species(from_info, reaction);
            register_reactant_species(to_info, reaction);
            register_product_species(product_species0);
            register_product_species(product_species1);

            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                world_->new_voxel_private(product_species0, from_coord));
            if (!new_mol0.second)
            {
                throw IllegalState("no place for the first product. 4");
            }
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world_->new_voxel_private(product_species1, to_coord));
            if (!new_mol1.second)
            {
                throw IllegalState("no place for the first product. 5");
            }

            reaction.products.push_back(
                reaction_type::particle_type(
                    new_mol0.first.first,
                    this->private_voxel2voxel(new_mol0.first.second)));
            reaction.products.push_back(
                reaction_type::particle_type(
                    new_mol1.first.first,
                    this->private_voxel2voxel(new_mol1.first.second)));
        }
    }
    else if(aserial == dloc)
    {
        if (bserial == cloc)
        {
            register_reactant_species(from_info, reaction);
            register_reactant_species(to_info, reaction);
            register_product_species(product_species0);
            register_product_species(product_species1);

            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                world_->new_voxel_private(product_species0, to_coord));
            if (!new_mol0.second)
            {
                throw IllegalState("no place for the first product. 6");
            }
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world_->new_voxel_private(product_species1, from_coord));
            if (!new_mol1.second)
            {
                throw IllegalState("no place for the first product. 7");
            }

            reaction.products.push_back(
                reaction_type::particle_type(
                    new_mol0.first.first,
                    this->private_voxel2voxel(new_mol0.first.second)));
            reaction.products.push_back(
                reaction_type::particle_type(
                    new_mol1.first.first,
                    this->private_voxel2voxel(new_mol1.first.second)));
        }
    }

    throw IllegalState("no place for the first product. 8");
}

/*
 * the First Order Reaction
 */
std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::apply_first_order_reaction_(
        const ReactionRule& reaction_rule, const LatticeWorld::particle_info info)
{
    const ReactionRule::product_container_type& products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;

    const Species& product_species0(*(products.begin()));
    switch(products.size()) {
        case 0:
            world_->remove_voxel_private(info.first);
            break;
        case 1:
            apply_a2b(info, product_species0, reaction);
            break;
        case 2:
            if (!apply_a2bc(info, product_species0,
                        (*(++products.begin())), reaction)) {
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

void LatticeSimulator::apply_a2b(
    const LatticeWorld::particle_info pinfo,
    const Species& product_species,
    reaction_type& reaction)
{
    const LatticeWorld::private_coordinate_type coord(pinfo.first);
    const std::string toloc(world_->get_molecule_info(product_species).loc);
    const std::string aserial(get_serial(coord));

    if (aserial == toloc)
    {
        register_reactant_species(pinfo, reaction);
        register_product_species(product_species);

        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel_private(product_species, coord));
        if (!new_mol.second)
        {
            throw IllegalState("no place for the product. 9");
        }
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol.first.first,
                this->private_voxel2voxel(new_mol.first.second)));
    }
    else
    {
        std::pair<LatticeWorld::private_coordinate_type, bool> neighbor(
                world_->check_neighbor_private(coord));
        const std::string nserial(get_serial(neighbor.first));
        if (nserial == toloc)
        {
            world_->remove_voxel_private(neighbor.first);
            register_reactant_species(pinfo, reaction);
            register_product_species(product_species);

            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world_->new_voxel_private(product_species, neighbor.first));
            if (!new_mol.second)
            {
                throw IllegalState("no place for the product. 10");
            }
            reaction.products.push_back(
                reaction_type::particle_type(
                    new_mol.first.first,
                    this->private_voxel2voxel(new_mol.first.second)));
        }
    }

}

bool LatticeSimulator::apply_a2bc(
    const LatticeWorld::particle_info pinfo,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    const LatticeWorld::private_coordinate_type coord(pinfo.first);
    const std::string toloc0(world_->get_molecule_info(product_species0).loc),
                      toloc1(world_->get_molecule_info(product_species1).loc);
    std::pair<LatticeWorld::private_coordinate_type, bool> neighbor(
            world_->check_neighbor_private(coord));
    const std::string aserial(get_serial(coord));
    const std::string nserial(get_serial(neighbor.first));

    if (aserial == toloc0)
    {
        if (nserial != toloc1)
        {
            return false;
            // throw IllegalState("no place for the product.");
            // TODO
        }
        register_reactant_species(pinfo, reaction);

        register_product_species(product_species0);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world_->new_voxel_private(product_species0, coord));
        if (!new_mol0.second)
        {
            throw IllegalState("no place for the first product. 12");
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
            throw IllegalState("no place for the second product. 13");
        }
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol1.first.first,
                this->private_voxel2voxel(new_mol1.first.second)));
        return true;
    }
    else if (nserial == toloc1)
    {
        if (nserial != toloc0)
        {
            return false;
            // throw IllegalState("no place for the product.");
            // TODO
        }
        register_reactant_species(pinfo, reaction);

        register_product_species(product_species0);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world_->new_voxel_private(product_species0, neighbor.first));
        if (!new_mol0.second)
        {
            throw IllegalState("no place for the first product. 15");
        }
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol0.first.first,
                this->private_voxel2voxel(new_mol0.first.second)));

        register_product_species(product_species1);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world_->new_voxel_private(product_species1, coord));
        if (!new_mol1.second)
        {
            throw IllegalState("no place for the second product. 16");
        }
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol1.first.first,
                this->private_voxel2voxel(new_mol1.first.second)));
        return true;
    }

    throw IllegalState("no place for the product.");
}

void LatticeSimulator::register_product_species(const Species& product_species)
{
    if (!world_->has_species(product_species))
    {
        new_species_.push_back(product_species);
    }
}

void LatticeSimulator::register_reactant_species(
        const LatticeWorld::particle_info pinfo, reaction_type reaction) const
{
    const MolecularTypeBase* mtype(world_->get_molecular_type_private(pinfo.first));
    const std::string location(
            mtype->location()->is_vacant() ? "" : mtype->location()->species().serial());
    reaction.reactants.push_back(
        reaction_type::particle_type(
            pinfo.second,
            Voxel(mtype->species(), world_->private2coord(pinfo.first),
                mtype->radius(), mtype->D(), location)));
    world_->remove_voxel_private(pinfo.first);
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
        const std::pair<LatticeWorld::private_coordinate_type, bool>
            neighbor(world_->move_to_neighbor(
                mtype, loc, (*mtype)[i], rng->uniform_int(0, 11)));

        if (!neighbor.second)
        {
            const LatticeWorld::particle_info info((*mtype)[i]);
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
