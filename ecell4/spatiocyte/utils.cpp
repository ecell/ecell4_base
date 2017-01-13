#include "utils.hpp"


namespace ecell4
{

namespace spatiocyte
{

const Real calculate_dimensional_factor(
    const VoxelPool* mt0, const VoxelPool* mt1,
    const boost::shared_ptr<const SpatiocyteWorld>& world)
{
    const Real voxel_radius(world->voxel_radius());
    const Real unit_area(world->unit_area());

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
        //     factor = 1. / (6 * sqrt(2.0) * Dtot * voxel_radius);
        // else
        //     factor = 1. / (6 * sqrt(2.0) * D_A * voxel_radius);
        factor = 1. / (6 * sqrt(2.0) * Dtot * voxel_radius);
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
        factor = sqrt(2.0) / (3 * D_A * voxel_radius);
        if (mt1->is_structure()) // B is Surface
        {
            factor *= unit_area;
        }
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::THREE)
    {
        factor = sqrt(2.0) / (3 * D_B * voxel_radius);
        if (mt0->is_structure()) // A is Surface
        {
            factor *= unit_area;
        }
    }
    else
        throw NotSupported("The dimension of a structure must be two or three.");
    return factor;
}

const Real calculate_alpha(const ReactionRule& rr, const boost::shared_ptr<SpatiocyteWorld>& world)
{
    const ReactionRule::reactant_container_type& reactants(rr.reactants());
    if (reactants.size() != 2)
        return 1.0;

    const Species species[2] = {reactants.at(0), reactants.at(1)};
    const MoleculeInfo info[2] = {
        world->get_molecule_info(species[0]),
        world->get_molecule_info(species[1])
    };
    VoxelPool* mt[2];
    bool is_created[2] = {false, false};
    for (int i(0); i < 2; ++i) {
        try
        {
            mt[i] = world->find_voxel_pool(species[i]);
        }
        catch(NotFound e)
        {
            VoxelPool* location(&(VacantType::getInstance()));
            if (info[i].loc != "") {
                try
                {
                    location = world->find_voxel_pool(Species(info[i].loc));
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
    const Real factor(calculate_dimensional_factor(mt[0], mt[1], boost::const_pointer_cast<const SpatiocyteWorld>(world)));
    for (int i(0); i < 2; ++i)
        if (is_created[i])
            delete mt[i];
    const Real alpha(1.0 / (factor * rr.k()));
    return alpha < 1.0 ? alpha : 1.0;
}

} // spatiocyte

} // ecell4
