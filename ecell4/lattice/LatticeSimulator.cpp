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


    const std::vector<ReactionRule>& rules(model_->reaction_rules());
    for (std::vector<ReactionRule>::const_iterator i(rules.begin());
        i != rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        if (rr.reactants().size() != 0)
        {
            continue;
        }
        const boost::shared_ptr<EventScheduler::Event>
            zeroth_order_reaction_event(
                create_zeroth_order_reaction_event(rr, world_->t()));
        scheduler_.add(zeroth_order_reaction_event);
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

void LatticeSimulator::register_events(const Species& sp)
{
    if (world_->find_molecular_type(sp)->with_voxels())
    {
        //TODO: Call steps only if sp is assigned not to StructureType.
        const boost::shared_ptr<EventScheduler::Event> step_event(
                create_step_event(sp, world_->t()));
        scheduler_.add(step_event);
    }

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

boost::shared_ptr<EventScheduler::Event> LatticeSimulator::create_step_event(
        const Species& species, const Real& t)
{
    boost::shared_ptr<EventScheduler::Event> event(
        new StepEvent(this, species, t, alpha_));
    return event;
}

boost::shared_ptr<EventScheduler::Event>
LatticeSimulator::create_zeroth_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<EventScheduler::Event> event(new ZerothOrderReactionEvent(
                this, reaction_rule, t));
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

/*
 * the Zeroth Order Reaction
 */
std::pair<bool, LatticeSimulator::reaction_type>
    LatticeSimulator::apply_zeroth_order_reaction_(
        const ReactionRule& reaction_rule)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;
    // return std::pair<bool, reaction_type>(false, reaction);

    for (ReactionRule::product_container_type::const_iterator
        i(reaction_rule.products().begin());
        i != reaction_rule.products().end(); ++i)
    {
        const Species& sp(*i);
        const LatticeWorld::molecule_info_type
            info(world_->get_molecule_info(sp));
        register_product_species(sp);

        while (true) //TODO: Avoid an inifinite loop
        {
            const LatticeWorld::coordinate_type
                coord(world_->rng()->uniform_int(0, world_->size() - 1));
            const Voxel v(
                sp, world_->coord2private(coord),
                info.radius, info.D, info.loc);

            if (world_->on_structure(v))
            {
                continue;
            }

            const std::pair<std::pair<ParticleID, Voxel>, bool>
                retval(world_->new_voxel_private(v));
            if (retval.second)
            {
                reaction.products.push_back(
                    reaction_type::particle_type(
                        retval.first.first,
                        this->private_voxel2voxel(retval.first.second)));
                break;
            }
        }
    }
    reactions_.push_back(reaction_rule);
    return std::pair<bool, reaction_type>(true, reaction);
}

std::pair<bool, LatticeSimulator::reaction_type> LatticeSimulator::attempt_reaction_(
    const LatticeWorld::particle_info_type info, LatticeWorld::coordinate_type to_coord,
    const Real& alpha)
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
        const Real P(k * factor * alpha);
        accp += P;
        if (accp > 1)
        {
            std::cerr << "The total acceptance probability [" << accp
                << "] exceeds 1 for '" << from_species.serial()
                << "' and '" << to_species.serial() << "'." << std::endl;
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

    const LatticeWorld::private_coordinate_type from_coord(
        world_->coord2private(p0.second.coordinate()));
    const LatticeWorld::private_coordinate_type to_coord(
        world_->coord2private(p1.second.coordinate()));

    const LatticeWorld::particle_info_type from_info(from_coord, p0.first);
    const LatticeWorld::particle_info_type to_info(to_coord, p1.first);

    switch (products.size())
    {
        case 1:
            apply_ab2c(from_info, to_info, *(products.begin()), reaction);
            break;
        case 2:
            apply_ab2cd(
                from_info, to_info,
                *(products.begin()), *(++(products.begin())), reaction);
            break;
        default:
            return std::pair<bool, reaction_type>(false, reaction);
    }

    reactions_.push_back(reaction_rule);
    return std::pair<bool, reaction_type>(true, reaction);
}

void LatticeSimulator::apply_ab2c(
    const LatticeWorld::particle_info_type from_info,
    const LatticeWorld::particle_info_type to_info,
    const Species& product_species,
    reaction_type& reaction)
{
    // A and B (from_info and to_info) become C (product_species)
    const std::string location(world_->get_molecule_info(product_species).loc);
    const std::string fserial(get_serial(from_info.first));
    const std::string floc(get_location(from_info.first));
    const std::string tserial(get_serial(to_info.first));
    const std::string tloc(get_location(to_info.first));

    if (tserial == location || tloc == location)
    {
        // B is on the location of C, or the location itself.
        // Place C at the coordinate of B, and remove A.
        register_reactant_species(from_info, reaction);
        register_reactant_species(to_info, reaction);
        register_product_species(product_species);

        if (tserial != location)
        {
            world_->remove_voxel_private(to_info.first);
        }

        world_->remove_voxel_private(from_info.first);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel_private(product_species, to_info.first));

        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol.first.first,
                this->private_voxel2voxel(new_mol.first.second)));
    }
    else if (fserial == location || floc == location)
    {
        // A is on the location of C, or the location itself.
        // Place C at the coordinate of A, and remove B.
        register_reactant_species(from_info, reaction);
        register_reactant_species(to_info, reaction);
        register_product_species(product_species);

        if (fserial != location)
        {
            world_->remove_voxel_private(from_info.first);
        }

        world_->remove_voxel_private(to_info.first);
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel_private(product_species, from_info.first));

        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol.first.first,
                this->private_voxel2voxel(new_mol.first.second)));
    }
    else
    {
        throw IllegalState("no place for the product.");
    }
}

// Not tested yet
void LatticeSimulator::apply_ab2cd(
    const LatticeWorld::particle_info_type from_info,
    const LatticeWorld::particle_info_type to_info,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    const LatticeWorld::private_coordinate_type from_coord(from_info.first);
    const LatticeWorld::private_coordinate_type to_coord(to_info.first);
    const std::string aserial(get_serial(from_coord));
    const std::string aloc(get_location(from_coord));
    const std::string bserial(get_serial(to_coord));
    const std::string bloc(get_location(to_coord));
    const std::string cloc(world_->get_molecule_info(product_species0).loc);
    const std::string dloc(world_->get_molecule_info(product_species1).loc);

    if (aserial == cloc || aloc == cloc)
    {
        if (bserial == dloc || bloc == dloc)
        {
            register_reactant_species(from_info, reaction);
            register_reactant_species(to_info, reaction);
            register_product_species(product_species0);
            register_product_species(product_species1);

            world_->remove_voxel_private(from_info.first);
            world_->remove_voxel_private(to_info.first);
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                world_->new_voxel_private(product_species0, from_coord));
            if (!new_mol0.second)
            {
                throw IllegalState("no place for the first product.");
            }
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world_->new_voxel_private(product_species1, to_coord));
            if (!new_mol1.second)
            {
                throw IllegalState("no place for the first product.");
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
    else if(aserial == dloc || aloc == dloc)
    {
        if (bserial == cloc || bloc == dloc)
        {
            register_reactant_species(from_info, reaction);
            register_reactant_species(to_info, reaction);
            register_product_species(product_species0);
            register_product_species(product_species1);

            world_->remove_voxel_private(from_info.first);
            world_->remove_voxel_private(to_info.first);
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                world_->new_voxel_private(product_species0, to_coord));
            if (!new_mol0.second)
            {
                throw IllegalState("no place for the first product.");
            }
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world_->new_voxel_private(product_species1, from_coord));
            if (!new_mol1.second)
            {
                throw IllegalState("no place for the first product.");
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

    throw IllegalState("no place for the first product.");
}

/*
 * the First Order Reaction
 */
std::pair<bool, LatticeSimulator::reaction_type>
    LatticeSimulator::apply_first_order_reaction_(
        const ReactionRule& reaction_rule,
        const LatticeSimulator::reaction_type::particle_type& p)
{
    const ReactionRule::product_container_type& products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;

    const LatticeWorld::private_coordinate_type coord(
        world_->coord2private(p.second.coordinate()));
    const LatticeWorld::particle_info_type info(coord, p.first);

    switch (products.size())
    {
        case 0:
            world_->remove_voxel_private(coord);
            break;
        case 1:
            apply_a2b(info, *(products.begin()), reaction);
            break;
        case 2:
            if (!apply_a2bc(
                info, *(products.begin()), (*(++products.begin())), reaction))
            {
                return std::pair<bool, reaction_type>(false, reaction);
            }
            break;
        default:
            return std::pair<bool, reaction_type>(false, reaction);
    }

    reactions_.push_back(reaction_rule);
    return std::pair<bool, reaction_type>(true, reaction);
}

void LatticeSimulator::apply_a2b(
    const LatticeWorld::particle_info_type pinfo,
    const Species& product_species,
    reaction_type& reaction)
{
    // A (pinfo) becomes B (product_species)
    const LatticeWorld::private_coordinate_type coord(pinfo.first);
    const std::string location(world_->get_molecule_info(product_species).loc);
    const std::string aserial(get_serial(coord));
    const std::string aloc(get_location(coord));

    if (aserial == location || aloc == location)
    {
        // A is on the location of B, or the location itself.
        register_reactant_species(pinfo, reaction);
        register_product_species(product_species);

        if (aserial != location)
        {
            world_->remove_voxel_private(pinfo.first);
        }

        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel_private(product_species, coord));

        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol.first.first,
                this->private_voxel2voxel(new_mol.first.second)));
    }
    else
    {
        // A is NOT on the location of B.
        // B must be released into a neighbor, which is the location of B
        std::pair<LatticeWorld::private_coordinate_type, bool>
            neighbor(world_->check_neighbor_private(coord));
        const std::string nserial(get_serial(neighbor.first));

        if (nserial == location)
        {
            // The neighbor is the location of B.
            // Place B at the neighbor, and remove A.
            register_reactant_species(pinfo, reaction);
            register_product_species(product_species);

            world_->remove_voxel_private(pinfo.first);
            world_->remove_voxel_private(neighbor.first);
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world_->new_voxel_private(product_species, neighbor.first));

            reaction.products.push_back(
                reaction_type::particle_type(
                    new_mol.first.first,
                    this->private_voxel2voxel(new_mol.first.second)));
        }
    }
}

bool LatticeSimulator::apply_a2bc(
    const LatticeWorld::particle_info_type pinfo,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    // A (pinfo) becomes B and C (product_species0 and product_species1)
    // At least, one of A and B must be placed at the neighbor.
    const LatticeWorld::private_coordinate_type coord(pinfo.first);
    const std::string
        bloc(world_->get_molecule_info(product_species0).loc),
        cloc(world_->get_molecule_info(product_species1).loc);
    std::pair<LatticeWorld::private_coordinate_type, bool>
        neighbor(world_->check_neighbor_private(coord));
    const std::string aserial(get_serial(coord));
    const std::string aloc(get_location(coord));
    const std::string nserial(get_serial(neighbor.first));
    const std::string nloc(get_location(neighbor.first));

    if (aserial == bloc || aloc == bloc)
    {
        // A is on the location of B, or the locaiton itself
        // C must be placed at the neighbor

        // if (nserial != cloc && nloc != cloc)
        if (nserial != cloc)
        {
            //TODO: C cannot be on the neighbor.
            return false;
            // throw IllegalState("no place for the product.");
        }

        register_reactant_species(pinfo, reaction);
        register_product_species(product_species0);
        register_product_species(product_species1);

        if (aserial != bloc)
        {
            world_->remove_voxel_private(pinfo.first);
        }
        // world_->remove_voxel_private(neighbor.first);

        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world_->new_voxel_private(product_species0, coord));
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world_->new_voxel_private(product_species1, neighbor.first));

        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol0.first.first,
                this->private_voxel2voxel(new_mol0.first.second)));
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol1.first.first,
                this->private_voxel2voxel(new_mol1.first.second)));
        return true;
    }
    else if (aserial == cloc || aloc == cloc)
    {
        // A is on the location of C, or the locaiton itself
        // B must be placed at the neighbor

        // if (nserial != bloc && nloc != bloc)
        if (nserial != bloc)
        {
            //TODO: B cannot be on the neighbor.
            return false;
            // throw IllegalState("no place for the product.");
        }

        register_reactant_species(pinfo, reaction);
        register_product_species(product_species0);
        register_product_species(product_species1);

        if (aserial != cloc)
        {
            world_->remove_voxel_private(pinfo.first);
        }
        world_->remove_voxel_private(neighbor.first);

        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world_->new_voxel_private(product_species0, neighbor.first));
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world_->new_voxel_private(product_species1, coord));

        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol0.first.first,
                this->private_voxel2voxel(new_mol0.first.second)));
        reaction.products.push_back(
            reaction_type::particle_type(
                new_mol1.first.first,
                this->private_voxel2voxel(new_mol1.first.second)));
        return true;
    }
    else
    {
        throw IllegalState("no place for the products.");
    }

    return false;
}

void LatticeSimulator::register_product_species(const Species& product_species)
{
    if (!world_->has_species(product_species))
    {
        new_species_.push_back(product_species);
    }
}

void LatticeSimulator::register_reactant_species(
        const LatticeWorld::particle_info_type pinfo, reaction_type reaction) const
{
    const MolecularTypeBase* mtype(world_->get_molecular_type_private(pinfo.first));
    const std::string location(
            mtype->location()->is_vacant() ? "" : mtype->location()->species().serial());
    reaction.reactants.push_back(
        reaction_type::particle_type(
            pinfo.second,
            Voxel(mtype->species(), world_->private2coord(pinfo.first),
                mtype->radius(), mtype->D(), location)));
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
    world_->set_t(time);
    top.second->fire(); // top.second->time_ is updated in fire()

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

    // std::cout << "[" << num_steps_ << "] sp=" << species.serial() << " : t=" << t()
    //     << " : alpha=" << alpha << " : walk." << std::endl;

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());

    MolecularTypeBase* mtype(world_->find_molecular_type(species));
    MolecularTypeBase* loc(mtype->location());

    Integer imax(rng->binomial(alpha, mtype->size()));
    if (imax != mtype->size())
    {
        mtype->shuffle(*rng);
    }
    //XXX: mtype->shuffle(*rng);

    if (!mtype->with_voxels())
    {
        throw NotSupported("MolecularType must be with voxels.");
    }

    Integer i(0);
    while (i < imax)
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
                retval(attempt_reaction_(info, to_coord, alpha));
            if (retval.first)
            {
                --i;
                --imax;

                if (imax > mtype->size())
                {
                    imax = mtype->size(); //XXX: for a dimerization
                }
            }
        }

        ++i;
    }
}

} // lattice

} // ecell4
