#include <iostream>

#include <ecell4/core/NetworkModel.hpp>
#include "GillespieSimulator.hpp"

using namespace ecell4;
using namespace ecell4::gillespie;


int main(int argc, char **argv)
{
    Species sp1("A"), sp2("B");
    ReactionRule rr1(create_unimolecular_reaction_rule(sp1, sp2, 5.001));

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
    world->save("test.h5");

    GillespieSimulator sim(model, world);

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
