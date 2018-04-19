#include "SpatiocyteReactions.hpp"
#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

// Utilities

const std::string get_serial(boost::shared_ptr<SpatiocyteWorld> world,
        const SpatiocyteWorld::coordinate_type coord)
{
    boost::shared_ptr<const VoxelPool> mtype(world->get_voxel_pool_at(coord));
    return mtype->is_vacant() ? "" : mtype->species().serial();
}

const std::string get_location(boost::shared_ptr<SpatiocyteWorld> world,
        const SpatiocyteWorld::coordinate_type coord)
{
    boost::shared_ptr<const VoxelPool> mtype(world->get_voxel_pool_at(coord));
    if (mtype->is_vacant())
        return "";
    boost::shared_ptr<const VoxelPool> ltype(mtype->location());
    return ltype->is_vacant() ? "" : ltype->species().serial();
}

static inline void
make_product(boost::shared_ptr<SpatiocyteWorld> world,
             ReactionInfo& rinfo,
             const Species& species,
             const SpatiocyteWorld::coordinate_type coord)
{
    if (world->has_species(species) && world->find_voxel_pool(species)->is_structure())
    {
        if (boost::optional<ParticleID> new_pid = world->new_voxel_structure(species, coord))
        {
            rinfo.add_product(ReactionInfo::Item(*new_pid, species, coord));
        }
    }
    else
    {
        if (boost::optional<ParticleID> new_pid = world->new_voxel(species, coord))
        {
            rinfo.add_product(ReactionInfo::Item(*new_pid, species, coord));
        }
    }
}

// Application of reactions

ReactionInfo
apply_a2b(boost::shared_ptr<SpatiocyteWorld> world,
          const ReactionInfo::Item& reactant_item,
          const Species& product_species)
{
    const Voxel::coordinate_type coord(reactant_item.voxel.coordinate);
    const std::string bloc(world->get_molecule_info(product_species).loc);
    const std::string aserial(get_serial(world, coord));
    const std::string aloc(get_location(world, coord));
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
            world->remove_voxel(coord);
        }

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            make_product(world, rinfo, product_species, coord);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            std::pair<ParticleID, Species> id_species_pair(world->get_voxel_at(coord));
            rinfo.add_product(ReactionInfo::Item(id_species_pair.first,
                                                 id_species_pair.second,
                                                 Voxel(coord)));
        }
    }
    else
    {
        // A is NOT on the location of B.
        // B must be released into a neighbor, which is the location of B
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(coord, bloc));

        if (neighbor.second)
        {
            // The neighbor is the location of B.
            // Place B at the neighbor, and remove A.
            rinfo.add_reactant(reactant_item);

            world->remove_voxel(coord);

            make_product(world, rinfo, product_species, neighbor.first);
        }
    }
    return rinfo;
}

ReactionInfo apply_a2bc(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::Item& reactant_item,
        const Species& product_species0,
        const Species& product_species1)
{
    // A (pinfo) becomes B and C (product_species0 and product_species1)
    // At least, one of A and B must be placed at the neighbor.
    const SpatiocyteWorld::coordinate_type coord(reactant_item.voxel.coordinate);
    const std::string bserial(product_species0.serial()),
                      cserial(product_species1.serial()),
                      bloc(world->get_molecule_info(product_species0).loc),
                      cloc(world->get_molecule_info(product_species1).loc);
    const std::string aserial(get_serial(world, coord));
    const std::string aloc(get_location(world, coord));

    ReactionInfo rinfo(world->t());

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // A is the locaiton of B,
        // or A is on the location of B,
        // or B is the location of A
        // C must be placed at the neighbor

        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(coord, cloc));
        const std::string nserial(get_serial(world, neighbor.first));
        const std::string nloc(get_location(world, neighbor.first));

        if (!neighbor.second)
        {
            //TODO: C cannot be on the neighbor.
            return rinfo;
        }

        rinfo.add_reactant(reactant_item);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of a new B-molecule
            world->remove_voxel(coord);
        }

        // No need to remove the neighbor because it's the location of C
        // world->remove_voxel(neighbor.first);

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            make_product(world, rinfo, product_species0, coord);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            std::pair<ParticleID, Species> id_species_pair(world->get_voxel_at(coord));
            rinfo.add_product(ReactionInfo::Item(id_species_pair.first,
                                                 id_species_pair.second,
                                                 Voxel(coord)));
        }

        // Place a new C-molecule at the neighbor
        make_product(world, rinfo, product_species1, neighbor.first);

        return rinfo;
    }
    else if (aserial == cloc || aloc == cloc || aloc == cserial)
    {
        // A is the locaiton of C,
        // or A is on the location of C,
        // or C is the location of A
        // B must be placed at the neighbor
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(coord, bloc));
        const std::string nserial(get_serial(world, neighbor.first));
        const std::string nloc(get_location(world, neighbor.first));

        if (!neighbor.second)
        {
            //TODO: B cannot be on the neighbor.
            return rinfo;
        }

        rinfo.add_reactant(reactant_item);

        if (aserial != cloc)
        {
            // Remove A once if A is not the location of a new C-molecule
            world->remove_voxel(coord);
        }

        // No need to remove the neighbor because it's the location of B
        // world->remove_voxel(neighbor.first);

        // Place a new B-molecule at the neighbor
        make_product(world, rinfo, product_species0, neighbor.first);

        if (aloc != cserial)
        {
            // Place a new C-molecule at the position of A
            make_product(world, rinfo, product_species1, coord);
        }
        else
        {
            // When C is the location of A, it's enough to remove A
            std::pair<ParticleID, Species> id_species_pair(world->get_voxel_at(coord));
            rinfo.add_product(ReactionInfo::Item(id_species_pair.first,
                                                 id_species_pair.second,
                                                 Voxel(coord)));
        }
        return rinfo;
    }
    return rinfo;
}

ReactionInfo apply_vanishment(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::Item& reactant_item0,
        const ReactionInfo::Item& reactant_item1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(reactant_item0);
    rinfo.add_reactant(reactant_item1);

    world->remove_voxel(reactant_item0.voxel.coordinate);
    world->remove_voxel(reactant_item1.voxel.coordinate);

    return rinfo;
}

ReactionInfo apply_ab2c(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::Item& reactant_item0,
        const ReactionInfo::Item& reactant_item1,
        const Species& product_species)
{
    const SpatiocyteWorld::coordinate_type coord0(reactant_item0.voxel.coordinate);
    const SpatiocyteWorld::coordinate_type coord1(reactant_item1.voxel.coordinate);

    // A and B (from_info and to_info) become C (product_species)
    const std::string location(world->get_molecule_info(product_species).loc);
    const std::string fserial(get_serial(world, coord0));
    const std::string floc(get_location(world, coord0));
    const std::string tserial(get_serial(world, coord1));
    const std::string tloc(get_location(world, coord1));

    ReactionInfo rinfo(world->t());

    if (tserial == location || tloc == location)
    {
        // B is on the location of C, or the location itself.
        // Place C at the coordinate of B, and remove A.
        rinfo.add_reactant(reactant_item0);
        rinfo.add_reactant(reactant_item1);

        if (tserial != location)
        {
            world->remove_voxel(coord1);
        }

        world->remove_voxel(coord0);

        make_product(world, rinfo, product_species, coord1);
    }
    else if (fserial == location || floc == location)
    {
        // A is on the location of C, or the location itself.
        // Place C at the coordinate of A, and remove B.
        rinfo.add_reactant(reactant_item0);
        rinfo.add_reactant(reactant_item1);

        if (fserial != location)
        {
            world->remove_voxel(coord0);
        }

        world->remove_voxel(coord1);

        make_product(world, rinfo, product_species, coord0);
    }
    return rinfo;
}

// For apply_ab2cd
ReactionInfo apply_ab2cd_in_order(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::Item& reactant_item0,
        const ReactionInfo::Item& reactant_item1,
        const Species& product_species0,
        const Species& product_species1,
        const SpatiocyteWorld::coordinate_type coord0,
        const SpatiocyteWorld::coordinate_type coord1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(reactant_item0);
    rinfo.add_reactant(reactant_item1);

    make_product(world, rinfo, product_species0, coord0);
    make_product(world, rinfo, product_species1, coord1);

    return rinfo;
}

ReactionInfo apply_ab2cd(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::Item& reactant_item0,
        const ReactionInfo::Item& reactant_item1,
        const Species& product_species0,
        const Species& product_species1)
{
    const SpatiocyteWorld::coordinate_type from_coord(reactant_item0.voxel.coordinate);
    const SpatiocyteWorld::coordinate_type to_coord(reactant_item1.voxel.coordinate);
    const std::string aserial(get_serial(world, from_coord));
    const std::string aloc(get_location(world, from_coord));
    const std::string bserial(get_serial(world, to_coord));
    const std::string bloc(get_location(world, to_coord));
    const std::string cloc(world->get_molecule_info(product_species0).loc);
    const std::string dloc(world->get_molecule_info(product_species1).loc);

    if (aserial == cloc || aloc == cloc)
    {
        if (bserial == dloc || bloc == dloc)
        {
            if (aserial != cloc)
            {
                // Remove A once if A is not the location of C
                world->remove_voxel(from_coord);
            }
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world->remove_voxel(to_coord);
            }
            return apply_ab2cd_in_order(
                world, reactant_item0, reactant_item1, product_species0, product_species1,
                from_coord, to_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world->check_neighbor(to_coord, dloc));

            if (neighbor.second)
            {
                world->remove_voxel(to_coord);
                if (aserial != cloc)
                {
                    // Remove A once if A is not the location of C
                    world->remove_voxel(from_coord);
                }
                return apply_ab2cd_in_order(
                    world, reactant_item0, reactant_item1, product_species0, product_species1,
                    from_coord, neighbor.first);
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
                world->remove_voxel(from_coord);
            }
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world->remove_voxel(to_coord);
            }
            return apply_ab2cd_in_order(
                world, reactant_item0, reactant_item1, product_species0, product_species1,
                to_coord, from_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world->check_neighbor(to_coord, cloc));

            if (neighbor.second)
            {
                world->remove_voxel(to_coord);
                if (aserial != dloc)
                {
                    // Remove A once if A is not the location of D
                    world->remove_voxel(from_coord);
                }
                return apply_ab2cd_in_order(
                    world, reactant_item0, reactant_item1, product_species0, product_species1,
                    neighbor.first, from_coord);
            }
        }
    }
    else if (bserial == cloc || bloc == cloc)
    {
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(to_coord, dloc));

        if (neighbor.second)
        {
            world->remove_voxel(from_coord);
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world->remove_voxel(to_coord);
            }
            return apply_ab2cd_in_order(
                world, reactant_item0, reactant_item1, product_species0, product_species1,
                to_coord, neighbor.first);
        }
    }
    else if (bserial == dloc || bloc == dloc)
    {
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(to_coord, dloc));

        if (neighbor.second)
        {
            world->remove_voxel(from_coord);
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world->remove_voxel(to_coord);
            }
            return apply_ab2cd_in_order(
                world, reactant_item0, reactant_item1, product_species0, product_species1,
                neighbor.first, to_coord);
        }
    }
    return ReactionInfo(world->t());
}

ReactionInfo apply_second_order_reaction(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionRule& reaction_rule,
        const ReactionInfo::Item& reactant_item0,
        const ReactionInfo::Item& reactant_item1)
{
    const ReactionRule::product_container_type& products(reaction_rule.products());

    switch (products.size())
    {
        case 0:
            return apply_vanishment(world, reactant_item0, reactant_item1);
        case 1:
            return apply_ab2c(world, reactant_item0, reactant_item1, *(products.begin()));
        case 2:
            return apply_ab2cd(world, reactant_item0, reactant_item1,
                               products.at(0), products.at(1));
        default:
            return ReactionInfo(world->t());
    }
}

} // spatiocyte

} // ecell4
