#include "SpatiocyteEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

/// ZerothOrderReactionEvent

ZerothOrderReactionEvent::ZerothOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), world_(world), rule_(rule)
{
    time_ = t + draw_dt();
}

void ZerothOrderReactionEvent::fire()
{
    reaction_info_type rinfo(world_->t());

    for (ReactionRule::product_container_type::const_iterator
        i(rule_.products().begin());
        i != rule_.products().end(); ++i)
    {
        const Species& sp(*i);
        const SpatiocyteWorld::molecule_info_type
            info(world_->get_molecule_info(sp));

        while (true) //TODO: Avoid an inifinite loop
        {
            // const SpatiocyteWorld::coordinate_type
            //     coord(world_->rng()->uniform_int(0, world_->size() - 1));
            const SpatiocyteWorld::coordinate_type
                coord(inner2coordinate(*world_,
                            world_->rng()->uniform_int(0, world_->inner_size() - 1)));
            const Voxel v(sp, coord, info.radius, info.D, info.loc);

            if (world_->on_structure(v))
            {
                continue;
            }

            const std::pair<std::pair<ParticleID, Voxel>, bool>
                retval(world_->new_voxel(v));
            if (retval.second)
            {
                rinfo.add_product(retval.first);
                break;
            }
        }
    }
    push_reaction(std::make_pair(rule_, rinfo));
    time_ += draw_dt();
}

Real ZerothOrderReactionEvent::draw_dt()
{
    const Real k(rule_.k());
    const Real p = k * world_->volume();
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(world_->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}


/// FirstOrderReactionEvent

FirstOrderReactionEvent::FirstOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), world_(world), rule_(rule)
{
    //assert(rule_.reactants().size() == 1);
    time_ = t + draw_dt();
}

void FirstOrderReactionEvent::fire()
{
    const reaction_info_type::particle_id_pair_type& p(
            world_->choice(*(rule_.reactants().begin())));
    const ReactionRule::product_container_type& products(rule_.products());

    time_ += draw_dt();
    switch (products.size())
    {
        case 0:
            {
                world_->remove_voxel(p.second.coordinate());
                reaction_info_type rinfo(world_->t());
                rinfo.add_reactant(p);
                push_reaction(std::make_pair(rule_, rinfo));
            }
            break;
        case 1:
            push_reaction(std::make_pair(rule_, apply_a2b(p, *(products.begin()))));
            break;
        case 2:
            {
                std::pair<bool, reaction_info_type>
                    retval(apply_a2bc(p, *(products.begin()), (*(++products.begin()))));
                if (retval.first)
                    push_reaction(std::make_pair(rule_, retval.second));
            }
            break;
    }
}

Real FirstOrderReactionEvent::draw_dt()
{
    const Species& reactant(*(rule_.reactants().begin()));
    const Integer num_r(world_->num_voxels_exact(reactant));
    const Real k(rule_.k());
    const Real p = k * num_r;
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(world_->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}

SpatiocyteEvent::reaction_info_type FirstOrderReactionEvent::apply_a2b(
    const ReactionInfo::particle_id_pair_type& p,
    const Species& product_species)
{
    // A (pinfo) becomes B (product_species)
    const SpatiocyteWorld::coordinate_type coord(p.second.coordinate());
    const std::string bloc(world_->get_molecule_info(product_species).loc);
    const std::string aserial(get_serial(coord));
    const std::string aloc(get_location(coord));
    const std::string bserial(product_species.serial());

    reaction_info_type rinfo(world_->t());

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // A is the location of B (B can be placed on A),
        // or A is on the location of B,
        // or A is on B.
        rinfo.add_reactant(p);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of B
            world_->remove_voxel(p.second.coordinate());
        }

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world_->new_voxel(product_species, coord));
            rinfo.add_product(new_mol.first);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            rinfo.add_product(world_->get_voxel(coord));
        }
    }
    else
    {
        // A is NOT on the location of B.
        // B must be released into a neighbor, which is the location of B
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world_->check_neighbor(coord, bloc));

        if (neighbor.second)
        {
            // The neighbor is the location of B.
            // Place B at the neighbor, and remove A.
            rinfo.add_reactant(p);

            world_->remove_voxel(p.second.coordinate());
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world_->new_voxel(product_species, neighbor.first));

            rinfo.add_product(new_mol.first);
        }
    }
    return rinfo;
}

std::pair<bool, SpatiocyteEvent::reaction_info_type>
FirstOrderReactionEvent::apply_a2bc(
    const SpatiocyteEvent::reaction_info_type::particle_id_pair_type& p,
    const Species& product_species0,
    const Species& product_species1)
{
    // A (pinfo) becomes B and C (product_species0 and product_species1)
    // At least, one of A and B must be placed at the neighbor.
    const SpatiocyteWorld::coordinate_type coord(p.second.coordinate());
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

        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world_->check_neighbor(coord, cloc));
        const std::string nserial(get_serial(neighbor.first));
        const std::string nloc(get_location(neighbor.first));

        if (!neighbor.second)
        {
            //TODO: C cannot be on the neighbor.
            return std::make_pair(false, reaction_info_type());
        }

        reaction_info_type rinfo(world_->t());
        rinfo.add_reactant(p);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of a new B-molecule
            world_->remove_voxel(p.second.coordinate());
        }

        // No need to remove the neighbor because it's the location of C
        // world_->remove_voxel(neighbor.first);

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                world_->new_voxel(product_species0, coord));
            rinfo.add_product(new_mol0.first);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            rinfo.add_product(world_->get_voxel(coord));
        }

        // Place a new C-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world_->new_voxel(product_species1, neighbor.first));
        rinfo.add_product(new_mol1.first);
        return std::make_pair(true, rinfo);
    }
    else if (aserial == cloc || aloc == cloc || aloc == cserial)
    {
        // A is the locaiton of C,
        // or A is on the location of C,
        // or C is the location of A
        // B must be placed at the neighbor
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world_->check_neighbor(coord, bloc));
        const std::string nserial(get_serial(neighbor.first));
        const std::string nloc(get_location(neighbor.first));

        if (!neighbor.second)
        {
            //TODO: B cannot be on the neighbor.
            return std::make_pair(false, reaction_info_type());
        }

        reaction_info_type rinfo(world_->t());
        rinfo.add_reactant(p);

        if (aserial != cloc)
        {
            // Remove A once if A is not the location of a new C-molecule
            world_->remove_voxel(p.second.coordinate());
        }

        // No need to remove the neighbor because it's the location of B
        // world_->remove_voxel(neighbor.first);

        // Place a new B-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world_->new_voxel(product_species0, neighbor.first));
        rinfo.add_product(new_mol0.first);

        if (aloc != cserial)
        {
            // Place a new C-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world_->new_voxel(product_species1, coord));
            rinfo.add_product(new_mol1.first);
        }
        else
        {
            // When C is the location of A, it's enough to remove A
            rinfo.add_product(world_->get_voxel(coord));
        }
        return std::make_pair(true, rinfo);
    }
    else
    {
        throw IllegalState(
            "no place for the products ["
            + product_species0.serial() + ", " + product_species1.serial() + "].");
    }
    return std::make_pair(false, reaction_info_type());
}

} // spatiocyte

} // ecell4
