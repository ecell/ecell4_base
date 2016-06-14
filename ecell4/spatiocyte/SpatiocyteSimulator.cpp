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
        const boost::shared_ptr<EventScheduler::Event>
            zeroth_order_reaction_event(
                create_zeroth_order_reaction_event(rr, world_->t()));
        scheduler_.add(zeroth_order_reaction_event);
    }

    nids_.clear();
    for (unsigned int i(0); i < 12; ++i)
        nids_.push_back(i);

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
            const Species sp(reactants.at(i));
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

boost::shared_ptr<EventScheduler::Event> SpatiocyteSimulator::create_step_event(
        const Species& species, const Real& t)
{
    double alpha(alpha_);
    alpha_map_type::const_iterator itr(alpha_map_.find(species));
    if (itr != alpha_map_.end() && (*itr).second < alpha)
        alpha = (*itr).second;

    boost::shared_ptr<EventScheduler::Event> event(
        new StepEvent(this, species, t, alpha));
    return event;
}

boost::shared_ptr<EventScheduler::Event>
SpatiocyteSimulator::create_zeroth_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<EventScheduler::Event> event(new ZerothOrderReactionEvent(
                this, reaction_rule, t));
    return event;
}

boost::shared_ptr<EventScheduler::Event>
SpatiocyteSimulator::create_first_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<EventScheduler::Event> event(new FirstOrderReactionEvent(
                this, reaction_rule, t));
    return event;
}

void SpatiocyteSimulator::finalize()
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
std::pair<bool, SpatiocyteSimulator::reaction_type>
    SpatiocyteSimulator::apply_zeroth_order_reaction_(
        const ReactionRule& reaction_rule)
{
    // const ReactionRule::product_container_type&
    //     products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;
    // return std::pair<bool, reaction_type>(false, reaction);

    for (ReactionRule::product_container_type::const_iterator
        i(reaction_rule.products().begin());
        i != reaction_rule.products().end(); ++i)
    {
        const Species& sp(*i);
        const SpatiocyteWorld::molecule_info_type
            info(world_->get_molecule_info(sp));
        register_product_species(sp);

        while (true) //TODO: Avoid an inifinite loop
        {
            const SpatiocyteWorld::coordinate_type
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
                    reaction_type::particle_type(retval.first));
                break;
            }
        }
    }
    last_reactions_.push_back(std::make_pair(reaction_rule, reaction_info_type(world_->t(), reaction.reactants, reaction.products)));
    return std::pair<bool, reaction_type>(true, reaction);
}

Real SpatiocyteSimulator::calculate_dimensional_factor(
    const MolecularTypeBase* mt0, const MolecularTypeBase* mt1) const
{
    const Species
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
        //     factor = 1. / (6 * sqrt(2.0) * Dtot * world_->voxel_radius());
        // else
        //     factor = 1. / (6 * sqrt(2.0) * D_A * world_->voxel_radius());
        factor = 1. / (6 * sqrt(2.0) * Dtot * world_->voxel_radius());
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
        factor = sqrt(2.0) / (3 * D_A * world_->voxel_radius());
        if (mt1->is_structure()) // B is Surface
        {
            factor *= world_->unit_area();
        }
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::THREE)
    {
        factor = sqrt(2.0) / (3 * D_B * world_->voxel_radius());
        if (mt0->is_structure()) // A is Surface
        {
            factor *= world_->unit_area();
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
    MolecularTypeBase* mt[2];
    bool is_created[2] = {false, false};
    for (int i(0); i < 2; ++i) {
        try
        {
            mt[i] = world_->find_molecular_type(species[i]);
        }
        catch(NotFound e)
        {
            MolecularTypeBase *location(&(VacantType::getInstance()));
            if (info[i].loc != "") {
                try
                {
                    location = world_->find_molecular_type(Species(info[i].loc));
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
    const Real factor(calculate_dimensional_factor(mt[0], mt[1]));
    for (int i(0); i < 2; ++i)
        if (is_created[i])
            delete mt[i];
    const Real alpha(1.0 / (factor * rule.k()));
    return alpha < 1.0 ? alpha : 1.0;
}

std::pair<SpatiocyteSimulator::attempt_reaction_result_type, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::attempt_reaction_(
    const SpatiocyteWorld::particle_info_type info, SpatiocyteWorld::coordinate_type to_coord,
    const Real& alpha)
{
    const MolecularTypeBase* from_mt(
        world_->get_molecular_type_private(info.first));
    const MolecularTypeBase* to_mt(
        world_->get_molecular_type_private(to_coord));

    if (to_mt->is_vacant())
    {
        return std::pair<attempt_reaction_result_type, reaction_type>(NO_REACTION, reaction_type());
    }

    const Species
        speciesA(from_mt->species()),
        speciesB(to_mt->species());

    const std::vector<ReactionRule> rules(
        model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return std::pair<attempt_reaction_result_type, reaction_type>(NO_REACTION, reaction_type());
    }

    const Real factor(calculate_dimensional_factor(from_mt, to_mt));

    const Real rnd(world_->rng()->uniform(0,1));
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
                << "] exceeds 1 for '" << speciesA.serial()
                << "' and '" << speciesB.serial() << "'." << std::endl;
        }
        if (accp >= rnd)
        {
            std::pair<bool, SpatiocyteSimulator::reaction_type>
                retval = apply_second_order_reaction_(
                    *itr,
                    world_->make_pid_voxel_pair(from_mt, info),
                    world_->make_pid_voxel_pair(to_mt, to_coord));
            if (retval.first)
            {
                return std::pair<attempt_reaction_result_type, reaction_type>(REACTION_SUCCEEDED, retval.second);
            }
            else
            {
                return std::pair<attempt_reaction_result_type, reaction_type>(REACTION_FAILED, reaction_type());
            }
        }
    }
    return std::pair<attempt_reaction_result_type, reaction_type>(REACTION_FAILED, reaction_type());
}

/*
 * the Reaction between two molecules
 */
std::pair<bool, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::apply_second_order_reaction_(
    const ReactionRule& reaction_rule,
    const SpatiocyteSimulator::reaction_type::particle_type& p0,
    const SpatiocyteSimulator::reaction_type::particle_type& p1)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;

    const SpatiocyteWorld::private_coordinate_type from_coord(
        world_->coord2private(p0.second.coordinate()));
    const SpatiocyteWorld::private_coordinate_type to_coord(
        world_->coord2private(p1.second.coordinate()));

    const SpatiocyteWorld::particle_info_type from_info(from_coord, p0.first);
    const SpatiocyteWorld::particle_info_type to_info(to_coord, p1.first);

    switch (products.size())
    {
        case 0:
            apply_vanishment(from_info, to_info, reaction);
            break;
        case 1:
            if (!apply_ab2c(from_info, to_info, *(products.begin()), reaction))
            {
                return std::pair<bool, reaction_type>(false, reaction);
            }
            break;
        case 2:
            if (!apply_ab2cd(
                    from_info, to_info,
                    *(products.begin()), *(++(products.begin())), reaction))
            {
                return std::pair<bool, reaction_type>(false, reaction);
            }
            break;
        default:
            return std::pair<bool, reaction_type>(false, reaction);
    }

    last_reactions_.push_back(std::make_pair(reaction_rule, reaction_info_type(world_->t(), reaction.reactants, reaction.products)));
    return std::pair<bool, reaction_type>(true, reaction);
}

void SpatiocyteSimulator::apply_vanishment(
    const SpatiocyteWorld::particle_info_type from_info,
    const SpatiocyteWorld::particle_info_type to_info,
    reaction_type& reaction)
{
    register_reactant_species(from_info, reaction);
    register_reactant_species(to_info, reaction);
    world_->remove_voxel_private(from_info.first);
    world_->remove_voxel_private(to_info.first);
}

bool SpatiocyteSimulator::apply_ab2c(
    const SpatiocyteWorld::particle_info_type from_info,
    const SpatiocyteWorld::particle_info_type to_info,
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
            reaction_type::particle_type(new_mol.first));
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
            reaction_type::particle_type(new_mol.first));
    }
    else
    {
        return false;
        // throw IllegalState(
        //     "no place for the product [" + product_species.serial() + "].");
    }
    return true;
}

// Not tested yet
bool SpatiocyteSimulator::apply_ab2cd(
    const SpatiocyteWorld::particle_info_type from_info,
    const SpatiocyteWorld::particle_info_type to_info,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    const SpatiocyteWorld::private_coordinate_type from_coord(from_info.first);
    const SpatiocyteWorld::private_coordinate_type to_coord(to_info.first);
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

            if (aserial != cloc)
            {
                // Remove A once if A is not the location of C
                world_->remove_voxel_private(from_info.first);
            }
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world_->remove_voxel_private(to_info.first);
            }
            apply_ab2cd_in_order(from_coord, product_species0,
                    to_coord, product_species1, reaction);
            return true;
        }
        else
        {
            std::pair<SpatiocyteWorld::private_coordinate_type, bool>
                neighbor(world_->check_neighbor_private(to_coord, dloc));

            if (neighbor.second)
            {
                register_reactant_species(from_info, reaction);
                register_reactant_species(to_info, reaction);
                register_product_species(product_species0);
                register_product_species(product_species1);

                world_->remove_voxel_private(to_info.first);
                if (aserial != cloc)
                {
                    // Remove A once if A is not the location of C
                    world_->remove_voxel_private(from_info.first);
                }
                apply_ab2cd_in_order(from_coord, product_species0,
                        neighbor.first, product_species1, reaction);
                return true;
            }
        }
    }
    else if (aserial == dloc || aloc == dloc)
    {
        if (bserial == cloc || bloc == dloc)
        {
            register_reactant_species(from_info, reaction);
            register_reactant_species(to_info, reaction);
            register_product_species(product_species0);
            register_product_species(product_species1);

            if (aserial != dloc)
            {
                // Remove A once if A is not the location of D
                world_->remove_voxel_private(from_info.first);
            }
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world_->remove_voxel_private(to_info.first);
            }
            apply_ab2cd_in_order(to_coord, product_species0,
                    from_coord, product_species1, reaction);
            return true;
        }
        else
        {
            std::pair<SpatiocyteWorld::private_coordinate_type, bool>
                neighbor(world_->check_neighbor_private(to_coord, cloc));

            if (neighbor.second)
            {
                register_reactant_species(from_info, reaction);
                register_reactant_species(to_info, reaction);
                register_product_species(product_species0);
                register_product_species(product_species1);

                world_->remove_voxel_private(to_info.first);
                if (aserial != dloc)
                {
                    // Remove A once if A is not the location of D
                    world_->remove_voxel_private(from_info.first);
                }
                apply_ab2cd_in_order(neighbor.first, product_species0,
                        from_coord, product_species1, reaction);
                return true;
            }
        }
    }
    else if (bserial == cloc || bloc == cloc)
    {
        std::pair<SpatiocyteWorld::private_coordinate_type, bool>
            neighbor(world_->check_neighbor_private(to_coord, dloc));

        if (neighbor.second)
        {
            register_reactant_species(from_info, reaction);
            register_reactant_species(to_info, reaction);
            register_product_species(product_species0);
            register_product_species(product_species1);

            world_->remove_voxel_private(from_info.first);
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world_->remove_voxel_private(to_info.first);
            }
            apply_ab2cd_in_order(to_coord, product_species0,
                    neighbor.first, product_species1, reaction);
            return true;
        }
    }
    else if (bserial == dloc || bloc == dloc)
    {
        std::pair<SpatiocyteWorld::private_coordinate_type, bool>
            neighbor(world_->check_neighbor_private(to_coord, dloc));

        if (neighbor.second)
        {
            register_reactant_species(from_info, reaction);
            register_reactant_species(to_info, reaction);
            register_product_species(product_species0);
            register_product_species(product_species1);

            world_->remove_voxel_private(from_info.first);
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world_->remove_voxel_private(to_info.first);
            }
            apply_ab2cd_in_order(neighbor.first, product_species0,
                    to_coord, product_species1, reaction);
            return true;
        }
    }
    // else
    // {
    //     throw IllegalState("Not Supported.");
    // }
    return false;
}

void SpatiocyteSimulator::apply_ab2cd_in_order(
    const SpatiocyteWorld::private_coordinate_type coord0,
    const Species& product_species0,
    const SpatiocyteWorld::private_coordinate_type coord1,
    const Species& product_species1,
    reaction_type& reaction)
{
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world_->new_voxel_private(product_species0, coord0));
        if (!new_mol0.second)
        {
            throw IllegalState("no place for " + product_species0.serial());
        }
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world_->new_voxel_private(product_species1, coord1));
        if (!new_mol1.second)
        {
            throw IllegalState("no place for " + product_species1.serial());
        }

        reaction.products.push_back(
            reaction_type::particle_type(new_mol0.first));
        reaction.products.push_back(
            reaction_type::particle_type(new_mol1.first));
}

/*
 * the First Order Reaction
 */
std::pair<bool, SpatiocyteSimulator::reaction_type>
    SpatiocyteSimulator::apply_first_order_reaction_(
        const ReactionRule& reaction_rule,
        const SpatiocyteSimulator::reaction_type::particle_type& p)
{
    const ReactionRule::product_container_type& products(reaction_rule.products());
    reaction_type reaction;
    reaction.rule = reaction_rule;

    const SpatiocyteWorld::private_coordinate_type coord(
        world_->coord2private(p.second.coordinate()));
    const SpatiocyteWorld::particle_info_type info(coord, p.first);

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

    last_reactions_.push_back(std::make_pair(reaction_rule, reaction_info_type(world_->t(), reaction.reactants, reaction.products)));
    return std::pair<bool, reaction_type>(true, reaction);
}

void SpatiocyteSimulator::apply_a2b(
    const SpatiocyteWorld::particle_info_type pinfo,
    const Species& product_species,
    reaction_type& reaction)
{
    // A (pinfo) becomes B (product_species)
    const SpatiocyteWorld::private_coordinate_type coord(pinfo.first);
    const std::string bloc(world_->get_molecule_info(product_species).loc);
    const std::string aserial(get_serial(coord));
    const std::string aloc(get_location(coord));
    const std::string bserial(product_species.serial());

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // A is the location of B (B can be placed on A),
        // or A is on the location of B,
        // or A is on B.
        register_reactant_species(pinfo, reaction);
        register_product_species(product_species);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of B
            world_->remove_voxel_private(pinfo.first);
        }

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world_->new_voxel_private(product_species, coord));
            reaction.products.push_back(
                reaction_type::particle_type(new_mol.first));
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            reaction.products.push_back(
                reaction_type::particle_type(
                    world_->get_voxel(world_->private2coord(coord))));
        }
    }
    else
    {
        // A is NOT on the location of B.
        // B must be released into a neighbor, which is the location of B
        std::pair<SpatiocyteWorld::private_coordinate_type, bool>
            neighbor(world_->check_neighbor_private(coord, bloc));

        if (neighbor.second)
        {
            // The neighbor is the location of B.
            // Place B at the neighbor, and remove A.
            register_reactant_species(pinfo, reaction);
            register_product_species(product_species);

            world_->remove_voxel_private(pinfo.first);
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world_->new_voxel_private(product_species, neighbor.first));

            reaction.products.push_back(
                reaction_type::particle_type(new_mol.first));
        }
    }
}

bool SpatiocyteSimulator::apply_a2bc(
    const SpatiocyteWorld::particle_info_type pinfo,
    const Species& product_species0,
    const Species& product_species1,
    reaction_type& reaction)
{
    // A (pinfo) becomes B and C (product_species0 and product_species1)
    // At least, one of A and B must be placed at the neighbor.
    const SpatiocyteWorld::private_coordinate_type coord(pinfo.first);
    const std::string
        bserial(product_species0.serial()),
        cserial(product_species1.serial()),
        bloc(world_->get_molecule_info(product_species0).loc),
        cloc(world_->get_molecule_info(product_species1).loc);
    const std::string aserial(get_serial(coord));
    const std::string aloc(get_location(coord));

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // A is the locaiton of B,
        // or A is on the location of B,
        // or B is the location of A
        // C must be placed at the neighbor

        std::pair<SpatiocyteWorld::private_coordinate_type, bool>
            neighbor(world_->check_neighbor_private(coord, cloc));
        const std::string nserial(get_serial(neighbor.first));
        const std::string nloc(get_location(neighbor.first));

        if (!neighbor.second)
        {
            //TODO: C cannot be on the neighbor.
            return false;
        }

        register_reactant_species(pinfo, reaction);
        register_product_species(product_species0);
        register_product_species(product_species1);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of a new B-molecule
            world_->remove_voxel_private(pinfo.first);
        }

        // No need to remove the neighbor because it's the location of C
        // world_->remove_voxel_private(neighbor.first);

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                world_->new_voxel_private(product_species0, coord));
            reaction.products.push_back(
                reaction_type::particle_type(new_mol0.first));
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            reaction.products.push_back(
                reaction_type::particle_type(
                    world_->get_voxel(world_->private2coord(coord))));
        }

        // Place a new C-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world_->new_voxel_private(product_species1, neighbor.first));
        reaction.products.push_back(
            reaction_type::particle_type(new_mol1.first));
        return true;
    }
    else if (aserial == cloc || aloc == cloc || aloc == cserial)
    {
        // A is the locaiton of C,
        // or A is on the location of C,
        // or C is the location of A
        // B must be placed at the neighbor
        std::pair<SpatiocyteWorld::private_coordinate_type, bool>
            neighbor(world_->check_neighbor_private(coord, bloc));
        const std::string nserial(get_serial(neighbor.first));
        const std::string nloc(get_location(neighbor.first));

        if (!neighbor.second)
        {
            //TODO: B cannot be on the neighbor.
            return false;
        }

        register_reactant_species(pinfo, reaction);
        register_product_species(product_species0);
        register_product_species(product_species1);

        if (aserial != cloc)
        {
            // Remove A once if A is not the location of a new C-molecule
            world_->remove_voxel_private(pinfo.first);
        }

        // No need to remove the neighbor because it's the location of B
        // world_->remove_voxel_private(neighbor.first);

        // Place a new B-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world_->new_voxel_private(product_species0, neighbor.first));
        reaction.products.push_back(
            reaction_type::particle_type(new_mol0.first));

        if (aloc != cserial)
        {
            // Place a new C-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world_->new_voxel_private(product_species1, coord));
            reaction.products.push_back(
                reaction_type::particle_type(new_mol1.first));
        }
        else
        {
            // When C is the location of A, it's enough to remove A
            reaction.products.push_back(
                reaction_type::particle_type(
                    world_->get_voxel(world_->private2coord(coord))));
        }
        return true;
    }
    else
    {
        throw IllegalState(
            "no place for the products ["
            + product_species0.serial() + ", " + product_species1.serial() + "].");
    }

    return false;
}

void SpatiocyteSimulator::register_product_species(const Species& product_species)
{
    if (!world_->has_species(product_species))
    {
        new_species_.push_back(product_species);
    }
}

void SpatiocyteSimulator::register_reactant_species(
        const SpatiocyteWorld::particle_info_type pinfo, reaction_type& reaction) const
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
    new_species_.clear();
    dt_ = scheduler_.next_time() - t();
    return false;
}

void SpatiocyteSimulator::step_()
{
    last_reactions_.clear();
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

    // update_alpha_map(); // may be performance cost
    for (std::vector<Species>::const_iterator itr(new_species_.begin());
        itr != new_species_.end(); ++itr)
    {
        register_events(*itr);
    }

    num_steps_++;
}

void SpatiocyteSimulator::walk(const Species& species)
{
    walk(species, 1.0);
}

void SpatiocyteSimulator::walk(const Species& species, const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    const MolecularTypeBase* mtype(world_->find_molecular_type(species));

    if (!mtype->with_voxels())
    {
        throw NotSupported("MolecularType must be with voxels.");
    }

    MolecularTypeBase::container_type voxels;
    copy(mtype->begin(), mtype->end(), back_inserter(voxels));

    if (mtype->get_dimension() == Shape::THREE)
        walk_in_space_(mtype, alpha);
    else // dimension == TWO, etc.
        walk_on_surface_(mtype, alpha);
}

void SpatiocyteSimulator::walk_in_space_(const MolecularTypeBase* mtype, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MolecularTypeBase::container_type voxels;
    copy(mtype->begin(), mtype->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (MolecularTypeBase::container_type::iterator itr(voxels.begin());
            itr != voxels.end(); ++itr)
    {
        const Integer rnd(rng->uniform_int(0, 11));
        const SpatiocyteWorld::particle_info_type info(*itr);
        if (world_->get_molecular_type_private(info.first) != mtype)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }
        const SpatiocyteWorld::private_coordinate_type neighbor(
                world_->get_neighbor_private_boundary(info.first, rnd));
        if (world_->can_move(info.first, neighbor))
        {
            if (rng->uniform(0,1) <= alpha)
                world_->move_private(info.first, neighbor, /*candidate=*/idx);
        }
        else
        {
            attempt_reaction_(info, neighbor, alpha);
        }
        ++idx;
    }
}

void SpatiocyteSimulator::walk_on_surface_(const MolecularTypeBase* mtype, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MolecularTypeBase::container_type voxels;
    copy(mtype->begin(), mtype->end(), back_inserter(voxels));

    const MolecularTypeBase* location(mtype->location());
    std::size_t idx(0);
    for (MolecularTypeBase::container_type::iterator itr(voxels.begin());
            itr != voxels.end(); ++itr)
    {
        const SpatiocyteWorld::particle_info_type info(*itr);
        if (world_->get_molecular_type_private(info.first) != mtype)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }

        ecell4::shuffle(*(rng.get()), nids_);
        for (std::vector<unsigned int>::const_iterator itr(nids_.begin());
                itr != nids_.end(); ++itr)
        {
            const SpatiocyteWorld::private_coordinate_type neighbor(
                    world_->get_neighbor_private_boundary(info.first, *itr));
            const MolecularTypeBase* target(world_->get_molecular_type_private(neighbor));

            if (target->get_dimension() > mtype->get_dimension())
            {
                continue;
            }

            if (world_->can_move(info.first, neighbor))
            {
                if (rng->uniform(0,1) <= alpha)
                    world_->move_private(info.first, neighbor, /*candidate=*/idx);
                break;
            }
            else
            {
                if (attempt_reaction_(info, neighbor, alpha).first != NO_REACTION)
                    break;
            }
        }
        ++idx;
    }
}

} // spatiocyte

} // ecell4
