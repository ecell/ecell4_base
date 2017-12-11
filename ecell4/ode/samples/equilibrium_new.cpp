#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/ode/ODESimulator.hpp>
//#include <ecell4/ode/ODENetworkModel.hpp>
#include <ecell4/ode/ODEReactionRule.hpp>
#include <ecell4/ode/ODESimulator_New.hpp>


using namespace ecell4;
using namespace ecell4::ode;

const bool use_coeff = true;

/**
 * main function
 */
int main(int argc, char** argv)
{
    const Real L(1e-6);
    const Real3 edge_lengths(L, L, L);
    const Real volume(L * L * L);
    const Real N(100);
    const Real ka(0.1), U(0.5);

    Species sp1("A"), sp2("B"), sp3("C");
    ReactionRule rr1, rr2;
    rr1.set_k(ka);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);
    if (use_coeff) {
        boost::shared_ptr<ReactionRuleDescriptorCPPfunc> rr1_desc(new ReactionRuleDescriptorCPPfunc(NULL) );
        std::vector<Real> rr1_left;
        std::vector<Real> rr1_right;
        rr1_left.push_back(2.0);
        rr1_right.push_back(2.0);
        rr1_right.push_back(2.0);
        rr1_desc->set_reactant_coefficients(rr1_left);
        rr1_desc->set_product_coefficients(rr1_right);
        rr1.set_descriptor(rr1_desc);
        if (rr1.has_descriptor() != true) {
            throw;
        }
    }

    const Real kd(ka * volume * (1 - U) / (U * U * N));
    rr2.set_k(kd);
    rr2.add_reactant(sp2);
    rr2.add_reactant(sp3);
    rr2.add_product(sp1);

    //if (use_coeff) {
    //    boost::shared_ptr<ReactionRuleDescriptorCPPfunc> rr2_desc(new ReactionRuleDescriptorCPPfunc(NULL) );
    //    std::vector<Real> rr2_left;
    //    std::vector<Real> rr2_right;
    //    rr2_left.push_back(2.0);
    //    rr2_left.push_back(2.0);
    //    rr2_right.push_back(2.0);
    //    rr2_desc->set_reactant_coefficients(rr2_left);
    //    rr2_desc->set_product_coefficients(rr2_right);
    //    rr2.set_descriptor(rr2_desc);
    //}

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);
    model->add_species_attribute(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    //boost::shared_ptr<ODENetworkModel>  ode_model(new ODENetworkModel(model) );
    ReactionRule rr3;
    rr3.add_reactant(sp1);
    rr3.add_product(sp3);

    //if (use_coeff) {
    //    boost::shared_ptr<ReactionRuleDescriptorCPPfunc> rr3_desc(new ReactionRuleDescriptorCPPfunc(NULL) );
    //    std::vector<Real> rr3_left;
    //    std::vector<Real> rr3_right;
    //    rr3_left.push_back(2.0);
    //    rr3_right.push_back(2.0);
    //    rr3_desc->set_reactant_coefficients(rr3_left);
    //    rr3_desc->set_product_coefficients(rr3_right);
    //    rr3.set_descriptor(rr3_desc);
    //}

    model->add_reaction_rule(rr3);
    //model->dump_reactions();

    boost::shared_ptr<ODEWorld_New> world(new ODEWorld_New(edge_lengths));
    world->add_molecules(sp1, N);

    ODESimulator_New target(model, world, RUNGE_KUTTA_CASH_KARP54 );
    target.initialize();

    Real next_time(0.0), dt(0.01);
    std::cout << target.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for (unsigned int i(0); i < 10000; ++i)
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
