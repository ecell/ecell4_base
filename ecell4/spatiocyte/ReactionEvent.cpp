#include "SpatiocyteEvent.hpp"
#include "SpatiocyteSimulator.hpp" // TODO: remove this include line

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
    SpatiocyteSimulator* sim, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), sim_(sim), world_(sim->world()), rule_(rule)
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
            {
                std::pair<bool, std::pair<ReactionRule, reaction_info_type> >
                    retval(sim_->apply_a2b(rule_, p, *(products.begin())));
                if (retval.first)
                    push_reaction(retval.second);
            }
            break;
        case 2:
            {
                std::pair<bool, std::pair<ReactionRule, reaction_info_type> >
                    retval(sim_->apply_a2bc(rule_, p, *(products.begin()), (*(++products.begin()))));
                if (retval.first)
                    push_reaction(retval.second);
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

} // spatiocyte

} // ecell4
