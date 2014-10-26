#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/lattice/LatticeSimulator.hpp>
typedef ecell4::lattice::LatticeWorld world_type;
typedef ecell4::lattice::LatticeSimulator simulator_type;

namespace ecell4
{

void run()
{
    const Real world_size(1);
    const Position3 edge_lengths(world_size, world_size, world_size);
    const Real voxel_radius(0.0025);

    const Integer N(60);

    const std::string D("1.0"), radius("0.0025");

    Species sp("A", radius, D);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(0);
    // rng->seed(time(NULL));

    boost::shared_ptr<world_type> world(
        new world_type(edge_lengths, voxel_radius, rng));

    std::cout << "col size = " << world->col_size() << ", row size = " << world->row_size() << ", layer size = " << world->layer_size() << std::endl;
    std::cout << "total size = " << world->size() << std::endl;

    world->add_molecules(sp, N);

    simulator_type sim(model, world);
    std::cout << "dt = " << sim.dt() << std::endl;
    while (sim.step(5.0));
}

} // ecell4

/**
 * main function
 */
int main(int argc, char** argv)
{
    ecell4::run();
}
