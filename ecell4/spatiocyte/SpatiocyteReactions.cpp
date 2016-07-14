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
    const VoxelPool* mtype(world->find_voxel_pool(coord));
    return mtype->is_vacant() ? "" : mtype->species().serial();
}

const std::string get_location(boost::shared_ptr<SpatiocyteWorld> world,
        const SpatiocyteWorld::coordinate_type coord)
{
    const VoxelPool* mtype(world->find_voxel_pool(coord));
    if (mtype->is_vacant())
        return "";
    const VoxelPool* ltype(mtype->location());
    return ltype->is_vacant() ? "" : ltype->species().serial();
}

// Application of reactions

ReactionInfo apply_a2b(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p,
        const Species& product_species)
{
    const SpatiocyteWorld::coordinate_type coord(p.second.coordinate());
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
        rinfo.add_reactant(p);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of B
            world->remove_voxel(p.second.coordinate());
        }

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world->new_voxel(product_species, coord));
            rinfo.add_product(new_mol.first);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            rinfo.add_product(world->get_voxel(coord));
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
            rinfo.add_reactant(p);

            world->remove_voxel(p.second.coordinate());
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world->new_voxel(product_species, neighbor.first));

            rinfo.add_product(new_mol.first);
        }
    }
    return rinfo;
}

ReactionInfo apply_a2bc(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p,
        const Species& product_species0,
        const Species& product_species1)
{
    // A (pinfo) becomes B and C (product_species0 and product_species1)
    // At least, one of A and B must be placed at the neighbor.
    const SpatiocyteWorld::coordinate_type coord(p.second.coordinate());
    const std::string
        bserial(product_species0.serial()),
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

        rinfo.add_reactant(p);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of a new B-molecule
            world->remove_voxel(p.second.coordinate());
        }

        // No need to remove the neighbor because it's the location of C
        // world->remove_voxel(neighbor.first);

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                    world->new_voxel(product_species0, coord));
            rinfo.add_product(new_mol0.first);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            rinfo.add_product(world->get_voxel(coord));
        }

        // Place a new C-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world->new_voxel(product_species1, neighbor.first));
        rinfo.add_product(new_mol1.first);
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

        rinfo.add_reactant(p);

        if (aserial != cloc)
        {
            // Remove A once if A is not the location of a new C-molecule
            world->remove_voxel(p.second.coordinate());
        }

        // No need to remove the neighbor because it's the location of B
        // world->remove_voxel(neighbor.first);

        // Place a new B-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world->new_voxel(product_species0, neighbor.first));
        rinfo.add_product(new_mol0.first);

        if (aloc != cserial)
        {
            // Place a new C-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world->new_voxel(product_species1, coord));
            rinfo.add_product(new_mol1.first);
        }
        else
        {
            // When C is the location of A, it's enough to remove A
            rinfo.add_product(world->get_voxel(coord));
        }
        return rinfo;
    }
    return rinfo;
}

ReactionInfo apply_vanishment(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(p0);
    rinfo.add_reactant(p1);

    world->remove_voxel(p0.second.coordinate());
    world->remove_voxel(p1.second.coordinate());

    return rinfo;
}

ReactionInfo apply_ab2c(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1,
        const Species& product_species)
{
    // A and B (from_info and to_info) become C (product_species)
    const std::string location(world->get_molecule_info(product_species).loc);
    const std::string fserial(get_serial(world, p0.second.coordinate()));
    const std::string floc(get_location(world, p0.second.coordinate()));
    const std::string tserial(get_serial(world, p1.second.coordinate()));
    const std::string tloc(get_location(world, p1.second.coordinate()));

    ReactionInfo rinfo(world->t());

    if (tserial == location || tloc == location)
    {
        // B is on the location of C, or the location itself.
        // Place C at the coordinate of B, and remove A.
        rinfo.add_reactant(p0);
        rinfo.add_reactant(p1);

        if (tserial != location)
        {
            world->remove_voxel(p1.second.coordinate());
        }

        world->remove_voxel(p0.second.coordinate());
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world->new_voxel(product_species, p1.second.coordinate()));

        rinfo.add_product(new_mol.first);
    }
    else if (fserial == location || floc == location)
    {
        // A is on the location of C, or the location itself.
        // Place C at the coordinate of A, and remove B.
        rinfo.add_reactant(p0);
        rinfo.add_reactant(p1);

        if (fserial != location)
        {
            world->remove_voxel(p0.second.coordinate());
        }

        world->remove_voxel(p1.second.coordinate());
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world->new_voxel(product_species, p0.second.coordinate()));

        rinfo.add_product(new_mol.first);
    }
    return rinfo;
}

// For apply_ab2cd
ReactionInfo apply_ab2cd_in_order(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1,
        const Species& product_species0,
        const Species& product_species1,
        const SpatiocyteWorld::coordinate_type coord0,
        const SpatiocyteWorld::coordinate_type coord1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(p0);
    rinfo.add_reactant(p1);

    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
        world->new_voxel(product_species0, coord0));
    if (!new_mol0.second)
    {
        throw IllegalState("no place for " + product_species0.serial());
    }
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
        world->new_voxel(product_species1, coord1));
    if (!new_mol1.second)
    {
        throw IllegalState("no place for " + product_species1.serial());
    }

    rinfo.add_product(new_mol0.first);
    rinfo.add_product(new_mol1.first);

    return rinfo;
}

ReactionInfo apply_ab2cd(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1,
        const Species& product_species0,
        const Species& product_species1)
{
    const SpatiocyteWorld::coordinate_type from_coord(p0.second.coordinate());
    const SpatiocyteWorld::coordinate_type to_coord(p1.second.coordinate());
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
                world->remove_voxel(p0.second.coordinate());
            }
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                from_coord, to_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world->check_neighbor(to_coord, dloc));

            if (neighbor.second)
            {
                world->remove_voxel(p1.second.coordinate());
                if (aserial != cloc)
                {
                    // Remove A once if A is not the location of C
                    world->remove_voxel(p0.second.coordinate());
                }
                return apply_ab2cd_in_order(
                    world, p0, p1, product_species0, product_species1,
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
                world->remove_voxel(p0.second.coordinate());
            }
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                to_coord, from_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world->check_neighbor(to_coord, cloc));

            if (neighbor.second)
            {
                world->remove_voxel(p1.second.coordinate());
                if (aserial != dloc)
                {
                    // Remove A once if A is not the location of D
                    world->remove_voxel(p0.second.coordinate());
                }
                return apply_ab2cd_in_order(
                    world, p0, p1, product_species0, product_species1,
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
            world->remove_voxel(p0.second.coordinate());
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                to_coord, neighbor.first);
        }
    }
    else if (bserial == dloc || bloc == dloc)
    {
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(to_coord, dloc));

        if (neighbor.second)
        {
            world->remove_voxel(p0.second.coordinate());
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                neighbor.first, to_coord);
        }
    }
    return ReactionInfo(world->t());
}

ReactionInfo apply_second_order_reaction(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionRule& reaction_rule,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());

    switch (products.size())
    {
        case 0:
            return apply_vanishment(world, p0, p1);
        case 1:
            return apply_ab2c(world, p0, p1, *(products.begin()));
        case 2:
            return apply_ab2cd(world, p0, p1,
                            *(products.begin()), *(++(products.begin())));
        default:
            return ReactionInfo(world->t());
    }
}

} // spatiocyte

} // ecell4
