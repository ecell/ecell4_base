#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/core/H5Save.hpp>

#include <../EGFRDSimulatorWrapper.hpp>
typedef ecell4::egfrd::EGFRDWorld world_type;
typedef ecell4::egfrd::EGFRDSimulatorWrapper simulator_type;

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

    const Real k2(ka), k1(kd);

    Species sp1("A", radius, D), sp2("B", radius, D), sp3("C", radius, D);
    ReactionRule rr1(create_unbinding_reaction_rule(sp1, sp2, sp3, k1)),
        rr2(create_binding_reaction_rule(sp2, sp3, sp1, k2));

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);
    model->add_species_attribute(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(time(NULL));

    const Integer matrix_size(3);
    boost::shared_ptr<world_type> world(
        new world_type(world_size, matrix_size, rng));

    // world->add_species(sp1);
    // world->add_species(sp2);
    // world->add_species(sp3);

    world->add_molecules(sp2, N);
    world->add_molecules(sp3, N);

    simulator_type sim(model, world);
    world->save("test_egfrd.h5");

    Real next_time(0.0), dt(0.02);
    std::cout << sim.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for (unsigned int i(0); i < N; ++i)
    {
        next_time += dt;
        while (sim.step(next_time)) {}

        std::cout << sim.t()
                  << "\t" << world->num_molecules(sp1)
                  << "\t" << world->num_molecules(sp2)
                  << "\t" << world->num_molecules(sp3)
                  << std::endl;
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
