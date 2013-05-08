#include <iostream>

#include <ecell4/core/NetworkModel.hpp>
#include "GillespieSimulator.hpp"

using namespace ecell4;
using namespace ecell4::gillespie;


int main(int argc, char **argv)
{
    Species sp1("A"), sp2("B");

    // Expand Reaction Rule.
    ReactionRule rr1;
    rr1.set_k(5.001);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_reaction_rule(rr1);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(time(NULL));

    Real volume(1.0);
    boost::shared_ptr<GillespieWorld> world(new GillespieWorld(volume, rng));
    // world->add_species(sp1);
    // world->add_species(sp2);
    world->add_molecules(sp1, 10);

    GillespieSimulator sim(model, world);
    sim.save_hdf5_init("test.h5");
    sim.save_hdf5();

    std::cout << "t = " << sim.t()
              << ", A: " << world->num_molecules(sp1)
              << ", B: " << world->num_molecules(sp2) << std::endl;
    for (int i = 0; i < 10; ++i)
    {
        sim.step();

        std::cout << "t = " << sim.t()
                  << ", A: " << world->num_molecules(sp1)
                  << ", B: " << world->num_molecules(sp2) << std::endl;
    }
}
