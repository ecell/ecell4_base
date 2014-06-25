#include <iostream>

#include <ecell4/core/types.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/gillespie/GillespieSimulator.hpp>

using namespace ecell4;
using namespace ecell4::gillespie;


int main(int argc, char **argv)
{
    Species sp1("A"), sp2("B"), sp3("C");
    const Real kf(0.25), kr(1.0);
    ReactionRule
        rr1(create_binding_reaction_rule(sp1, sp2, sp3, kf)),
        rr2(create_unbinding_reaction_rule(sp3, sp1, sp2, kr));

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);
    model->add_species_attribute(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    rng->seed(time(NULL));

    const Real L(1.0);
    const Position3 edge_lengths(L, L, L);
    boost::shared_ptr<GillespieWorld> world(new GillespieWorld(edge_lengths, rng));
    world->add_molecules(sp3, 10);
    world->save("test_gillespie.h5");

    GillespieSimulator sim(model, world);

    std::cout << "t = " << sim.t()
              << ", A: " << world->num_molecules(sp1)
              << ", B: " << world->num_molecules(sp2)
              << ", C: " << world->num_molecules(sp3) << std::endl;
    for (int i = 0; i < 100; ++i)
    {
        sim.step();

        std::cout << "t = " << sim.t()
                  << ", A: " << world->num_molecules(sp1)
                  << ", B: " << world->num_molecules(sp2)
                  << ", C: " << world->num_molecules(sp3) << std::endl;
    }

    return 0;
}
