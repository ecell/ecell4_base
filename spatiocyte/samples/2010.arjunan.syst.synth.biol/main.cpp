#include <string>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../../utils.hpp"
#include "../../SpatiocyteSimulator.hpp"

#include "../../CoordinateLogger.hpp"
#include "../../SpatiocyteVisualizationLogger.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;


const Real L = 1e-6;
const Real volume = L * L * L;
const Position3 edge_lengths = Position3(L, L, L);
const Real voxel_radius = 2.5e-9;
const Integer N = 60;

void setup_model(ecell4::NetworkModel& model)
{
    ecell4::Species Surface("Surface");

    ecell4::Species
        MinDatp("MinD(p=atp,bs,loc=cyt)", "16e-12"),
        MinDadp("MinD(p=adp,bs,loc=cyt)", "16e-12"),
        MinEE("MinEE(bs1,bs2,loc=cyt)", "10e-12"),
        Surface_MinD("MinD(p=atp,bs,loc=mem)", "0.02e-12"),
        Surface_MinEE("MinEE(bs1,bs2,loc=mem)", "0.02e-12"),
        Surface_MinDEE("MinD(p=atp,bs[1],loc=mem).MinEE(bs1[1],bs2,loc=mem)", "0.02e-12"),
        Surface_MinDEED("MinD(p=atp,bs[1],loc=mem).MinEE(bs1[1],bs2[2],loc=mem).MinD(p=atp,bs[2],loc=mem)", "0.02e-12");

    model.add_species(MinDatp);
    model.add_species(MinDadp);
    model.add_species(MinEE);
    model.add_species(Surface_MinD);
    model.add_species(Surface_MinEE);
    model.add_species(Surface_MinDEE);
    model.add_species(Surface_MinDEED);

    // model.add_reaction_rule(
    //     create_binding_reaction_rule(
    //         Surface, MinDatp, SurfaceMinD, 2.2e-8)); // rr1
    // model.add_reaction_rule(
    //     create_bibi_reaction_rule(
    //         Surface_MinD, MinDatp, MinD, Mind, 3e-20)); // rr2

    model.add_reaction_rule(
        create_binding_reaction_rule(
            Surface_MinD, MinEE, Surface_MinDEE, 5e-19)); // rr3
    model.add_reaction_rule(
        create_unbinding_reaction_rule(
            Surface_MinDEE, Surface_MinEE, MinDadp, 1)); // rr4
    model.add_reaction_rule(
        create_unimolecular_reaction_rule(MinDadp, MinDatp, 5)); // rr5
    model.add_reaction_rule(
        create_binding_reaction_rule(
            Surface_MinD, Surface_MinDEE, Surface_MinDEED, 5e-15)); // rr6
    model.add_reaction_rule(
        create_unbinding_reaction_rule(
            Surface_MinDEED, Surface_MinDEE, MinDadp, 1)); // rr7
    model.add_reaction_rule(
        create_unimolecular_reaction_rule(Surface_MinEE, MinEE, 0.83)); // rr8
}

/**
 * main function
 */
int main(int argc, char** argv)
{
    boost::shared_ptr<ecell4::NetworkModel> model(new ecell4::NetworkModel());
    setup_model(*model);

    boost::shared_ptr<SpatiocyteWorld> world(
        new SpatiocyteWorld(edge_lengths, voxel_radius));
    world->add_molecules((*model).species("MinD(p=adp,bs,loc=cyt)"), 1300);
    world->add_molecules((*model).species("MinD(p=atp,bs[1],loc=mem).MinEE(bs1[1],bs2,loc=mem)"), 700);

    SpatiocyteSimulator sim(model, world);

    SpatiocyteVisualizationLogger logger(world);

    logger.add_species((*model).species("MinEE(bs1,bs2,loc=mem)"));
    logger.add_species((*model).species("MinD(p=atp,bs[1],loc=mem).MinEE(bs1[1],bs2,loc=mem)"));
    logger.add_species((*model).species("MinD(p=atp,bs[1],loc=mem).MinEE(bs1[1],bs2[2],loc=mem).MinD(p=atp,bs[2],loc=mem)"));
    logger.add_species((*model).species("MinD(p=atp,bs,loc=mem)"));

    // Real next_time(0.0), dt(0.5);
    logger.initialize();
    logger.log();

    sim.save_hdf5_init( std::string("spatiocyte.hdf5") );

    Real next_time(0.0), dt(0.02);

    for (unsigned int i(0); i < 100; ++i)
    {
        next_time += dt;
        while (sim.step(next_time)) {}

        sim.save_hdf5();

        std::cout << sim.t()
                  << "\t" << world->num_molecules((*model).species("MinEE(bs1,bs2,loc=mem)"))
                  << "\t" << world->num_molecules((*model).species("MinD(p=atp,bs[1],loc=mem).MinEE(bs1[1],bs2,loc=mem)"))
                  << "\t" << world->num_molecules((*model).species("MinD(p=atp,bs[1],loc=mem).MinEE(bs1[1],bs2[2],loc=mem).MinD(p=atp,bs[2],loc=mem)"))
                  << std::endl;
    }

}
