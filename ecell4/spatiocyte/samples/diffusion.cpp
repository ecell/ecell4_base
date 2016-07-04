#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/spatiocyte/SpatiocyteWorld.hpp>
#include <ecell4/spatiocyte/SpatiocyteSimulator.hpp>
typedef ecell4::spatiocyte::SpatiocyteWorld world_type;
typedef ecell4::spatiocyte::SpatiocyteSimulator simulator_type;

namespace ecell4
{

void run()
{
    const Real world_size(1);
    const Real3 edge_lengths(world_size, world_size, world_size);
    const Real voxel_radius(0.0025);

    const Integer N(60);

    const std::string D("1.0"), radius("0.0025");

    Species sp("A", radius, D);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<RandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(0);
    // rng->seed(time(NULL));

    // boost::shared_ptr<world_type> world(
    //     new world_type(edge_lengths, voxel_radius, rng));
    // boost::shared_ptr<world_type> world(
    //     create_spatiocyte_world_vector_impl(edge_lengths, voxel_radius, rng));
    boost::shared_ptr<world_type> world(
        ecell4::spatiocyte::create_spatiocyte_world_cell_list_impl(
            edge_lengths, voxel_radius, Integer3(5, 5, 5), rng));

    std::cout << "col size = " << world->col_size()
        << ", row size = " << world->row_size()
        << ", layer size = " << world->layer_size() << std::endl;
    std::cout << "total size = " << world->inner_size() << std::endl;

    world->add_molecules(sp, N);

    simulator_type sim(model, world);
    std::cout << "dt = " << sim.dt() << std::endl;
    for (unsigned int i(0); i != 1000; ++i)
    {
        sim.step();
    }

    // while (sim.step(1.0)) ; // do nothing
}

} // ecell4

/**
 * main function
 */
int main(int argc, char** argv)
{
    ecell4::run();
}
