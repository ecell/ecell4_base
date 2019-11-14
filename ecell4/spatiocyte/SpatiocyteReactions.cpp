#include "SpatiocyteReactions.hpp"
#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

// Utilities

inline const std::string get_serial(boost::shared_ptr<SpatiocyteWorld> world,
                                    const Voxel &voxel)
{
    boost::shared_ptr<const VoxelPool> mtype(voxel.get_voxel_pool());
    return mtype->is_vacant() ? "" : mtype->species().serial();
}

inline const std::string get_location(boost::shared_ptr<SpatiocyteWorld> world,
                                      const Voxel &voxel)
{
    boost::shared_ptr<const VoxelPool> mtype(voxel.get_voxel_pool());
    if (mtype->is_vacant())
        return "";
    boost::shared_ptr<const VoxelPool> ltype(mtype->location());
    return ltype->is_vacant() ? "" : ltype->species().serial();
}

static inline void make_product(boost::shared_ptr<SpatiocyteWorld> world,
                                ReactionInfo &rinfo, const Species &species,
                                const Voxel voxel)
{
    if (world->has_species(species) &&
        world->find_voxel_pool(species)->is_structure())
    {
        if (boost::optional<ParticleID> new_pid =
                world->new_voxel_structure(species, voxel))
        {
            rinfo.add_product(ReactionInfo::Item(*new_pid, species, voxel));
        }
    }
    else
    {
        if (boost::optional<ParticleID> new_pid =
                world->new_particle(species, voxel))
        {
            rinfo.add_product(ReactionInfo::Item(*new_pid, species, voxel));
        }
    }
}

// Application of reactions

ReactionInfo apply_a2b(boost::shared_ptr<SpatiocyteWorld> world,
                       const ReactionInfo::Item &reactant_item,
                       const Species &product_species)
{
    const Voxel voxel(reactant_item.voxel);
    const std::string bloc(world->get_molecule_info(product_species).loc);
    const std::string aserial(get_serial(world, voxel));
    const std::string aloc(get_location(world, voxel));
    const std::string bserial(product_species.serial());

    ReactionInfo rinfo(world->t());

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // A is the location of B (B can be placed on A),
        // or A is on the location of B,
        // or A is on B.
        rinfo.add_reactant(reactant_item);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of B
            voxel.clear();
        }

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            make_product(world, rinfo, product_species, voxel);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            std::pair<ParticleID, Species> id_species_pair(
                world->get_voxel_at(voxel));
            rinfo.add_product(ReactionInfo::Item(
                id_species_pair.first, id_species_pair.second, voxel));
        }
    }
    else
    {
        // A is NOT on the location of B.
        // B must be released into a neighbor, which is the location of B
        if (boost::optional<Voxel> neighbor =
                world->check_neighbor(voxel, bloc))
        {
            // The neighbor is the location of B.
            // Place B at the neighbor, and remove A.
            rinfo.add_reactant(reactant_item);

            voxel.clear();

            make_product(world, rinfo, product_species, *neighbor);
        }
    }
    return rinfo;
}

ReactionInfo apply_a2bc(boost::shared_ptr<SpatiocyteWorld> world,
                        const ReactionInfo::Item &reactant_item,
                        const Species &product_species0,
                        const Species &product_species1)
{
    // A (pinfo) becomes B and C (product_species0 and product_species1)
    // At least, one of A and B must be placed at the neighbor.
    const Voxel voxel(reactant_item.voxel);
    const std::string bserial(product_species0.serial()),
        cserial(product_species1.serial()),
        bloc(world->get_molecule_info(product_species0).loc),
        cloc(world->get_molecule_info(product_species1).loc);
    const std::string aserial(get_serial(world, voxel));
    const std::string aloc(get_location(world, voxel));

    ReactionInfo rinfo(world->t());

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // A is the locaiton of B,
        // or A is on the location of B,
        // or B is the location of A
        // C must be placed at the neighbor

        boost::optional<Voxel> neighbor(world->check_neighbor(voxel, cloc));

        if (!neighbor)
        {
            // TODO: C cannot be on the neighbor.
            return rinfo;
        }

        const std::string nserial(get_serial(world, *neighbor));
        const std::string nloc(get_location(world, *neighbor));

        rinfo.add_reactant(reactant_item);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of a new B-molecule
            voxel.clear();
        }

        // No need to remove the neighbor because it's the location of C
        // neighbor.first.clear();

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            make_product(world, rinfo, product_species0, voxel);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            std::pair<ParticleID, Species> id_species_pair(
                world->get_voxel_at(voxel));
            rinfo.add_product(ReactionInfo::Item(
                id_species_pair.first, id_species_pair.second, voxel));
        }

        // Place a new C-molecule at the neighbor
        make_product(world, rinfo, product_species1, *neighbor);

        return rinfo;
    }
    else if (aserial == cloc || aloc == cloc || aloc == cserial)
    {
        // A is the locaiton of C,
        // or A is on the location of C,
        // or C is the location of A
        // B must be placed at the neighbor
        boost::optional<Voxel> neighbor(world->check_neighbor(voxel, bloc));

        if (!neighbor)
        {
            // TODO: B cannot be on the neighbor.
            return rinfo;
        }

        const std::string nserial(get_serial(world, *neighbor));
        const std::string nloc(get_location(world, *neighbor));

        rinfo.add_reactant(reactant_item);

        if (aserial != cloc)
        {
            // Remove A once if A is not the location of a new C-molecule
            voxel.clear();
        }

        // No need to remove the neighbor because it's the location of B
        // neighbor.first.clear();

        // Place a new B-molecule at the neighbor
        make_product(world, rinfo, product_species0, *neighbor);

        if (aloc != cserial)
        {
            // Place a new C-molecule at the position of A
            make_product(world, rinfo, product_species1, voxel);
        }
        else
        {
            // When C is the location of A, it's enough to remove A
            std::pair<ParticleID, Species> id_species_pair(
                world->get_voxel_at(voxel));
            rinfo.add_product(ReactionInfo::Item(
                id_species_pair.first, id_species_pair.second, voxel));
        }
        return rinfo;
    }
    return rinfo;
}

ReactionInfo apply_vanishment(boost::shared_ptr<SpatiocyteWorld> world,
                              const ReactionInfo::Item &reactant_item0,
                              const ReactionInfo::Item &reactant_item1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(reactant_item0);
    rinfo.add_reactant(reactant_item1);

    reactant_item0.voxel.clear();
    reactant_item1.voxel.clear();

    return rinfo;
}

ReactionInfo apply_ab2c(boost::shared_ptr<SpatiocyteWorld> world,
                        const ReactionInfo::Item &reactant_item0,
                        const ReactionInfo::Item &reactant_item1,
                        const Species &product_species)
{
    const Voxel voxel0(reactant_item0.voxel);
    const Voxel voxel1(reactant_item1.voxel);

    // A and B (from_info and to_info) become C (product_species)
    const std::string location(world->get_molecule_info(product_species).loc);
    const std::string fserial(get_serial(world, voxel0));
    const std::string floc(get_location(world, voxel0));
    const std::string tserial(get_serial(world, voxel1));
    const std::string tloc(get_location(world, voxel1));

    ReactionInfo rinfo(world->t());

    if (tserial == location || tloc == location)
    {
        // B is on the location of C, or the location itself.
        // Place C at the coordinate of B, and remove A.
        rinfo.add_reactant(reactant_item0);
        rinfo.add_reactant(reactant_item1);

        if (tserial != location)
        {
            voxel1.clear();
        }

        voxel0.clear();

        make_product(world, rinfo, product_species, voxel1);
    }
    else if (fserial == location || floc == location)
    {
        // A is on the location of C, or the location itself.
        // Place C at the coordinate of A, and remove B.
        rinfo.add_reactant(reactant_item0);
        rinfo.add_reactant(reactant_item1);

        if (fserial != location)
        {
            voxel0.clear();
        }

        voxel1.clear();

        make_product(world, rinfo, product_species, voxel0);
    }
    return rinfo;
}

// For apply_ab2cd
ReactionInfo apply_ab2cd_in_order(boost::shared_ptr<SpatiocyteWorld> world,
                                  const ReactionInfo::Item &reactant_item0,
                                  const ReactionInfo::Item &reactant_item1,
                                  const Species &product_species0,
                                  const Species &product_species1,
                                  const Voxel &voxel0, const Voxel &voxel1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(reactant_item0);
    rinfo.add_reactant(reactant_item1);

    make_product(world, rinfo, product_species0, voxel0);
    make_product(world, rinfo, product_species1, voxel1);

    return rinfo;
}

ReactionInfo apply_ab2cd(boost::shared_ptr<SpatiocyteWorld> world,
                         const ReactionInfo::Item &reactant_item0,
                         const ReactionInfo::Item &reactant_item1,
                         const Species &product_species0,
                         const Species &product_species1)
{
    const Voxel &src(reactant_item0.voxel);
    const Voxel &dst(reactant_item1.voxel);

    const std::string aserial(get_serial(world, src));
    const std::string aloc(get_location(world, src));
    const std::string bserial(get_serial(world, dst));
    const std::string bloc(get_location(world, dst));
    const std::string cloc(world->get_molecule_info(product_species0).loc);
    const std::string dloc(world->get_molecule_info(product_species1).loc);

    if (aserial == cloc || aloc == cloc)
    {
        if (bserial == dloc || bloc == dloc)
        {
            if (aserial != cloc)
            {
                // Remove A once if A is not the location of C
                src.clear();
            }
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                dst.clear();
            }
            return apply_ab2cd_in_order(world, reactant_item0, reactant_item1,
                                        product_species0, product_species1, src,
                                        dst);
        }
        else
        {
            boost::optional<Voxel> neighbor(world->check_neighbor(dst, dloc));

            if (neighbor)
            {
                dst.clear();
                if (aserial != cloc)
                {
                    // Remove A once if A is not the location of C
                    src.clear();
                }
                return apply_ab2cd_in_order(world, reactant_item0,
                                            reactant_item1, product_species0,
                                            product_species1, src, *neighbor);
            }
        }
    }
    else if (aserial == dloc || aloc == dloc)
    {
        if (bserial == cloc || bloc == dloc)
        {
            if (aserial != dloc)
            {
                // Remove A once if A is not the location of D
                src.clear();
            }
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                dst.clear();
            }
            return apply_ab2cd_in_order(world, reactant_item0, reactant_item1,
                                        product_species0, product_species1, dst,
                                        src);
        }
        else
        {
            boost::optional<Voxel> neighbor(world->check_neighbor(dst, cloc));

            if (neighbor)
            {
                dst.clear();
                if (aserial != dloc)
                {
                    // Remove A once if A is not the location of D
                    src.clear();
                }
                return apply_ab2cd_in_order(world, reactant_item0,
                                            reactant_item1, product_species0,
                                            product_species1, *neighbor, src);
            }
        }
    }
    else if (bserial == cloc || bloc == cloc)
    {
        if (boost::optional<Voxel> neighbor =
                (world->check_neighbor(dst, dloc)))
        {
            src.clear();
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                dst.clear();
            }
            return apply_ab2cd_in_order(world, reactant_item0, reactant_item1,
                                        product_species0, product_species1, dst,
                                        *neighbor);
        }
    }
    else if (bserial == dloc || bloc == dloc)
    {
        if (boost::optional<Voxel> neighbor = world->check_neighbor(dst, dloc))
        {
            src.clear();
            if (bserial != dloc)
            {
                // remove b once if b is not the location of d
                dst.clear();
            }
            return apply_ab2cd_in_order(world, reactant_item0, reactant_item1,
                                        product_species0, product_species1,
                                        *neighbor, dst);
        }
    }
    return ReactionInfo(world->t());
}

ReactionInfo
apply_second_order_reaction(boost::shared_ptr<SpatiocyteWorld> world,
                            const ReactionRule &reaction_rule,
                            const ReactionInfo::Item &reactant_item0,
                            const ReactionInfo::Item &reactant_item1)
{
    const ReactionRule::product_container_type &products(
        reaction_rule.products());

    switch (products.size())
    {
    case 0:
        return apply_vanishment(world, reactant_item0, reactant_item1);
    case 1:
        return apply_ab2c(world, reactant_item0, reactant_item1,
                          *(products.begin()));
    case 2:
        return apply_ab2cd(world, reactant_item0, reactant_item1,
                           products.at(0), products.at(1));
    default:
        return ReactionInfo(world->t());
    }
}

} // namespace spatiocyte

} // namespace ecell4
