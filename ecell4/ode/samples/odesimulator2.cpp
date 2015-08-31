
#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/ode/ODENetworkModel.hpp>
#include <ecell4/ode/ODEReactionRule.hpp>
#include <ecell4/ode/ODERatelaw.hpp>
#include <ecell4/ode/ODESimulator.hpp>
#include <boost/format.hpp>


using namespace ecell4;
using namespace ecell4::ode;

double f(ODERatelaw::state_container_type const &reac, ODERatelaw::state_container_type const &prod, 
        double const v, double const time, ODEReactionRule const &rr)
{
    return 0.0;
}

int main(int argc, char **argv)
{
    const Real L(1e-6);
    const Real3 edge_lengths(L, L, L);
    const Real volume(L * L * L);
    const Real N(60);
    const Real ka(0.1), U(0.5);

    Species sp1("A"), sp2("B"), sp3("C");
    ODEReactionRule rr1, rr2;
    rr1.add_reactant(sp1, 1.0);
    rr1.add_product(sp2, 1.0);
    rr1.add_product(sp3, 1.0);
    boost::shared_ptr<ODERatelawMassAction> ratelaw1(new ODERatelawMassAction(ka));
    boost::shared_ptr<ODERatelawCppCallback> ratelawcb(new ODERatelawCppCallback(f) );
    rr1.set_ratelaw(ratelaw1);
    std::cout << rr1.k() << std::endl;

    const Real kd(ka * volume * (1 - U) / (U * U * N));
    rr2.add_reactant(sp2, 1.0);
    rr2.add_reactant(sp3, 1.0);
    rr2.add_product(sp1, 1.0);
    boost::shared_ptr<ODERatelawMassAction> ratelaw2(new ODERatelawMassAction(kd));
    rr2.set_ratelaw(ratelaw2);

    boost::shared_ptr<ODENetworkModel> model(new ODENetworkModel());
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);
    model->dump_reactions();

    boost::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));
    world->add_molecules(sp1, N);

    ODESimulator sim(model, world, Controlled_Rosenbrock);
    sim.initialize();
    Real next_time(0.0), dt(0.01);
    std::cout << sim.t() 
              << "\t" << world->num_molecules(sp1) 
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for(unsigned int i(0); i < 100000; i++)
    {
        next_time += dt;
        sim.step(next_time);
        std::cout << sim.t() 
                  << "\t" << world->num_molecules(sp1) 
                  << "\t" << world->num_molecules(sp2)
                  << "\t" << world->num_molecules(sp3)
                  << std::endl;
    }
    return 0;
}
