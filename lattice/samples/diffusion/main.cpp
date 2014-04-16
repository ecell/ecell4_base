#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <../LatticeSimulator.hpp>
typedef ecell4::lattice::LatticeWorld world_type;
typedef ecell4::lattice::LatticeSimulator simulator_type;

namespace ecell4
{

void run()
{
    const Real world_size(1e-6);
    const Position3 edge_lengths(world_size, world_size, world_size);
    const Real volume(world_size * world_size * world_size);
    const Real voxel_radius(2.5e-9);

    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    Species sp("A", radius, D);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(time(NULL));

    boost::shared_ptr<world_type> world(
        new world_type(edge_lengths, voxel_radius, rng));

    world->add_molecules(sp, N);

    simulator_type sim(model, world);
    while(sim.step(1.0));
}

} // ecell4

/**
 * main function
 */
int main(int argc, char** argv)
{
    ecell4::run();
}
