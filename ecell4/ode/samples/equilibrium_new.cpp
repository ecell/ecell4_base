#include <iostream>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/ode/ODESimulator.hpp>


using namespace ecell4;
using namespace ecell4::ode;

const bool use_coeff = true;

// Real mass_action(const ReactionRuleDescriptor::state_container_type& r, const ReactionRuleDescriptor::state_container_type& p, Real volume, Real t, const ReactionRuleDescriptorCPPfunc& rd)
// {
//     Real ret = volume;
//     // Real ret = k_ * volume;
//     ReactionRuleDescriptor::state_container_type::const_iterator i(r.begin());
//     ReactionRuleDescriptor::reaction_coefficient_list_type::const_iterator j(rd.reactant_coefficients().begin());
//     for (; i != r.end() && j != rd.reactant_coefficients().end(); ++i, ++j)
//     {
//         ret *= std::pow((*i) / volume, (*j));
//     }
//     return ret;
// }

/**
 * main function
 */
int main(int argc, char** argv)
{
    const Real L(1e-6);
    const Real3 edge_lengths(L, L, L);
    const Real volume(L * L * L);
    const Real N(100);
    const Real U(0.5);
    const Real ka(use_coeff == 0 ? 0.1 : 0.1 / (N * U / volume));
    const Real kd(use_coeff == 0 ? ka * volume * (1 - U) / (U * U * N) : ka * std::pow((2 * U / (1 - U)), 2));

    Species sp1("A"), sp2("B"), sp3("C");

    ReactionRule rr1;
    rr1.set_k(ka);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);

    if (use_coeff)
    {
        // std::shared_ptr<ReactionRuleDescriptor> rr1_desc(new ReactionRuleDescriptorCPPfunc(mass_action));
        std::shared_ptr<ReactionRuleDescriptor> rr1_desc(new ReactionRuleDescriptorMassAction(rr1.k()));
        std::vector<Real> rr1_left;
        std::vector<Real> rr1_right;
        rr1_left.push_back(2.0);
        rr1_right.push_back(1.0);
        rr1_right.push_back(1.0);
        rr1_desc->set_reactant_coefficients(rr1_left);
        rr1_desc->set_product_coefficients(rr1_right);
        rr1_desc->resize_reactants(1);
        rr1_desc->resize_products(2);
        rr1.set_descriptor(rr1_desc);

        assert(rr1.has_descriptor());
    }

    ReactionRule rr2;
    rr2.set_k(kd);
    rr2.add_reactant(sp2);
    rr2.add_reactant(sp3);
    rr2.add_product(sp1);

    if (use_coeff)
    {
        // std::shared_ptr<ReactionRuleDescriptor> rr2_desc(new ReactionRuleDescriptorCPPfunc(mass_action));
        std::shared_ptr<ReactionRuleDescriptor> rr2_desc(new ReactionRuleDescriptorMassAction(rr2.k()));
        std::vector<Real> rr2_left;
        std::vector<Real> rr2_right;
        rr2_left.push_back(1.0);
        rr2_left.push_back(1.0);
        rr2_right.push_back(2.0);
        rr2_desc->set_reactant_coefficients(rr2_left);
        rr2_desc->set_product_coefficients(rr2_right);
        rr2_desc->resize_reactants(2);
        rr2_desc->resize_products(1);
        rr2.set_descriptor(rr2_desc);
    }

    std::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    // ReactionRule rr3;
    // rr3.add_reactant(sp1);
    // rr3.add_product(sp3);

    //if (use_coeff) {
    //    std::shared_ptr<ReactionRuleDescriptorCPPfunc> rr3_desc(new ReactionRuleDescriptorCPPfunc(NULL) );
    //    std::vector<Real> rr3_left;
    //    std::vector<Real> rr3_right;
    //    rr3_left.push_back(2.0);
    //    rr3_right.push_back(2.0);
    //    rr3_desc->set_reactant_coefficients(rr3_left);
    //    rr3_desc->set_product_coefficients(rr3_right);
    //    rr3.set_descriptor(rr3_desc);
    //}

    // model->add_reaction_rule(rr3);

    // model->dump_reactions();

    std::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));
    world->add_molecules(sp1, N);

    ODESimulator target(world, model, RUNGE_KUTTA_CASH_KARP54);
    target.initialize();

    std::cout << world->evaluate(model->reaction_rules()[0]) << std::endl;
    std::cout << world->evaluate(model->reaction_rules()[1]) << std::endl;

    Real next_time(0.0), dt(1.0);
    std::cout << target.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for (unsigned int i(0); i < 100; ++i)
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
