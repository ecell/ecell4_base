/**
   $ PREFIX=/foo/bar LD_LIBRARY_PATH=${PREFIX}/lib:${SRCPATH}/epdp \
   LIBRARY_PATH=${PREFIX}/lib:${SRCPATH}/epdp \
   CPLUS_INCLUDE_PATH=${PREFIX}/include:${SRCPATH}/epdp \
   g++ -DHAVE_CONFIG_H main.cpp -lgsl -lgslcblas -lgfrd \
   -lecell4-core -lecell4-ode -lecell4-gillespie -lecell4-egfrd -lecell4-bd
   $ PREFIX=/foo/bar LD_LIBRARY_PATH=${PREFIX}/lib:${SRCPATH}/epdp ./a.out
 */

#define EGFRD_MODE 0
#define BD_MODE 1
#define ODE_MODE 2
#define GILLESPIE_MODE 3

#ifndef STYPE
#define STYPE EGFRD_MODE
#endif

#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/core/H5Save.hpp>

#if STYPE == EGFRD_MODE
// #include <ecell4/egfrd/EGFRDSimulatorWrapper.hpp>
#include <../EGFRDSimulatorWrapper.hpp>
typedef ecell4::egfrd::EGFRDWorld world_type;
typedef ecell4::egfrd::EGFRDSimulatorWrapper simulator_type;
#elif STYPE == BD_MODE
#include <ecell4/bd/BDSimulator.hpp>
typedef ecell4::bd::BDWorld world_type;
typedef ecell4::bd::BDSimulator simulator_type;
#elif STYPE == ODE_MODE
#include <ecell4/ode/ODESimulator.hpp>
typedef ecell4::ode::ODEWorld world_type;
typedef ecell4::ode::ODESimulator simulator_type;
#else // STYPE == GILLESPIE_MODE
#include <ecell4/gillespie/GillespieSimulator.hpp>
typedef ecell4::gillespie::GillespieWorld world_type;
typedef ecell4::gillespie::GillespieSimulator simulator_type;
#endif

#include <hdf5.h>
#include <H5Cpp.h>

using namespace H5;

namespace ecell4
{

void run()
{
    const Real world_size(1e-6);
    const Position3 edge_lengths(world_size, world_size, world_size);
    const Real volume(world_size * world_size * world_size);

    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");
    // const Real kD(
    //     4 * M_PI * (2 * std::atof(D.c_str())) * (2 * std::atof(radius.c_str())));

    const Real kd(0.1), U(0.5);
    const Real ka(kd * volume * (1 - U) / (U * U * N));

#if (STYPE == EGFRD_MODE) || (STYPE == BD_MODE)
    const Real k2(ka), k1(kd);
#else
    const Real k2(ka * kD / (ka + kD));
    const Real k1(k2 * kd / ka);
#endif

    Species sp1("A", radius, D), sp2("B", radius, D), sp3("C", radius, D);
    ReactionRule rr1(create_unbinding_reaction_rule(sp1, sp2, sp3, k1)),
        rr2(create_binding_reaction_rule(sp2, sp3, sp1, k2));

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(time(NULL));

#if STYPE == EGFRD_MODE
    const Integer matrix_size(3);
    boost::shared_ptr<world_type> world(
        new world_type(world_size, matrix_size, rng));
#elif STYPE == BD_MODE
    boost::shared_ptr<world_type> world(new world_type(edge_lengths, rng));
#elif STYPE == ODE_MODE
    boost::shared_ptr<world_type> world(new world_type(volume));
#else // STYPE == GILLESPIE_MODE
    boost::shared_ptr<world_type> world(new world_type(volume, rng));
#endif

    world->add_species(sp1);
    world->add_species(sp2);
    world->add_species(sp3);

    world->add_molecules(sp1, N);

    simulator_type sim(model, world);
    // ecell4_hdf5_manager<world_type, int>
    //     hdf("hoge.h5", model, world, "ParticleSpace");

#if STYPE == BD_MODE
    sim.set_dt(1e-3);
#endif

    Real next_time(0.0), dt(0.02);
	sim.save_hdf5_init(std::string("mapk.hdf5"));
    std::cout << sim.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for (unsigned int i(0); i < 100; ++i)
    {
        next_time += dt;
        while (sim.step(next_time)) {}

        std::cout << sim.t()
                  << "\t" << world->num_molecules(sp1)
                  << "\t" << world->num_molecules(sp2)
                  << "\t" << world->num_molecules(sp3)
                  << std::endl;
		sim.save_hdf5();
        // hdf.save();
    }
}

} // ecell4

/**
 * main function
 */
int main(int argc, char** argv)
{
    ecell4::run();
}
