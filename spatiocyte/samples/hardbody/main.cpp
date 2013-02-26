#include <string>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../../utils.hpp"
#include "../../SpatiocyteSimulator.hpp"
#include "../../CoordinateLogger.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;

/**
 * main function
 */
int main(int argc, char** argv)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real volume(L * L * L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    const Real k2(0.1), U(0.5);
    const Real k1(k2 * volume * (1 - U) / (U * U * N));
    ecell4::Species sp1("A", radius, D), sp2("B", radius, D), sp3("C", radius, D);
    ReactionRule rr1(create_binding_reaction_rule(sp1, sp2, sp3, k1)),
        rr2(create_unbinding_reaction_rule(sp3, sp1, sp2, k2));

    boost::shared_ptr<ecell4::NetworkModel> model(new ecell4::NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<SpatiocyteWorld> world(
        new SpatiocyteWorld(edge_lengths, voxel_radius));
    world->add_species(sp1);
    world->add_species(sp2);
    world->add_species(sp3);
    world->add_molecules(sp3, N);

    SpatiocyteSimulator sim(model, world);
    sim.initialize();

    // assert(world->volume(), actual_volume_cuboid(edge_lengths, voxel_radius));

    CoordinateLogger logger(world);
    logger.add_species(sp3);

    logger.initialize();
    logger.log();
    std::cout << sim.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;

    Real next_time(0.0), dt(0.02);
    for (unsigned int i(0); i < 1000; ++i)
    {
        next_time += dt;
        while (sim.step(next_time)) {}

        logger.log();
        std::cout << sim.t()
                  << "\t" << world->num_molecules(sp1)
                  << "\t" << world->num_molecules(sp2)
                  << "\t" << world->num_molecules(sp3)
                  << std::endl;
    }
}
