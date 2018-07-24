
#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
// #include <ecell4/ode/ODENetworkModel.hpp>
// #include <ecell4/ode/ODEReactionRule.hpp>
// #include <ecell4/ode/ODERatelaw.hpp>
#include <ecell4/ode/ODESimulator.hpp>
#include <boost/format.hpp>


using namespace ecell4;
using namespace ecell4::ode;

double f(ReactionRuleDescriptor::state_container_type const &reac, ReactionRuleDescriptor::state_container_type const &prod, double const v, double const time, ReactionRuleDescriptorCPPfunc const &rr)
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
    ReactionRule rr1, rr2;
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);
    boost::shared_ptr<ReactionRuleDescriptorMassAction> ratelaw1(new ReactionRuleDescriptorMassAction(ka));
    ratelaw1->resize_reactants(1);
    ratelaw1->resize_products(2);
    rr1.set_descriptor(ratelaw1);
    std::cout << ratelaw1->k() << std::endl;

    const Real kd(ka * volume * (1 - U) / (U * U * N));
    rr2.add_reactant(sp2);
    rr2.add_reactant(sp3);
    rr2.add_product(sp1);
    boost::shared_ptr<ReactionRuleDescriptorMassAction> ratelaw2(new ReactionRuleDescriptorMassAction(kd));
    ratelaw2->resize_reactants(2);
    ratelaw2->resize_products(1);
    rr2.set_descriptor(ratelaw2);
    // std::cout << ratelaw2->k() << std::endl;

    boost::shared_ptr<Model> model(new NetworkModel());
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    for (Model::reaction_rule_container_type::const_iterator i(model->reaction_rules().begin());
        i != model->reaction_rules().end(); ++i)
    {
        std::cout << (*i).as_string() << std::endl;
    }

    boost::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));
    world->add_molecules(sp1, N);

    ODESimulator sim(model, world, ROSENBROCK4_CONTROLLER);
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
