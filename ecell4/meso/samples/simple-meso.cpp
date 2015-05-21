#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/core/Sphere.hpp>

#include <ecell4/meso/MesoscopicWorld.hpp>
#include <ecell4/meso/MesoscopicSimulator.hpp>
typedef ecell4::meso::MesoscopicWorld world_type;
typedef ecell4::meso::MesoscopicSimulator simulator_type;


namespace ecell4
{

void run()
{
    const Real L(10);
    const Real L_2(L * 0.5);
    const Real3 edge_lengths(L, L, L);
    const Integer3 matrix_sizes(30, 30, 30);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(Species("A", "0.0025", "1", "C"));

    boost::shared_ptr<RandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(0);
    // rng->seed(time(NULL));

    boost::shared_ptr<world_type> world(
        new world_type(edge_lengths, matrix_sizes, rng));
    world->bind_to(model);

    world->add_structure(
        Species("C"),
        boost::shared_ptr<const Shape>(
            new Sphere(Real3(L_2, L_2, L_2), L_2 * 0.8)));

    world->add_molecules(Species("A"), 1800);

    simulator_type sim(model, world);
    sim.run(1.0);
}

} // ecell4

/**
 * main function
 */
int main(int argc, char** argv)
{
    ecell4::run();
}
