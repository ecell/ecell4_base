#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/ode/ODESimulator_ratelow.hpp>

using namespace ecell4;
using namespace ecell4::ode;

double calculate_flux(double *state_array, double volume)
{
    double k = 1.0;
    double flux = k / volume;
    for(int i(0); i < 1; i++) {
        flux *= double(*(state_array + i)) * volume;
    }
    return flux;
}

/**
 * main function
 */
int main(int argc, char** argv)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);

    Species sp1("A"), sp2("B"), sp3("C");
    ReactionRule rr1;
    rr1.set_k(1.0);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);
    boost::shared_ptr<Ratelow> ratelow(new RatelowCppCallback(calculate_flux) );
    rr1.set_ratelow(ratelow);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);
    model->add_species_attribute(sp3);
    model->add_reaction_rule(rr1);

    boost::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));
    world->add_molecules(sp1, 60);
    world->save("test_ode.h5");

    world->bind_to(model);

    ODESimulatorRatelow target(model, world);

    Real next_time(0.0), dt(0.01);

    std::cout << target.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for (unsigned int i(0); i < 200; ++i)
    {
        next_time += dt;
        target.step(next_time);
        std::cout << target.t()
                  << "\t" << world->num_molecules(sp1)
                  << "\t" << world->num_molecules(sp2)
                  << "\t" << world->num_molecules(sp3)
                  << std::endl;
    }
}
