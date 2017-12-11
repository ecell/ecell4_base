
#include "ODESimulator_New.hpp"

#include <boost/numeric/odeint.hpp>
#include <algorithm>
#include <functional>

namespace odeint = boost::numeric::odeint;

namespace ecell4
{

namespace ode
{


void ODEWorld_New::bind_to(boost::shared_ptr<Model> model)
{
    //if (generated_)
    //{
    //    std::cerr << "Warning: NetworkModel is already bound to ODEWorld."
    //        << std::endl;
    //}
    //else if (model_.expired())
    //{
    //    std::cerr << "Warning: ODENetworkModel is already bound to ODEWorld."
    //        << std::endl;
    //}
    //try
    //{
    //    boost::shared_ptr<NetworkModel> tmp(new NetworkModel(model));
    //    generated_.swap(tmp);
    //    model_.reset();
    //}
    //catch (NotSupported e)
    //{
    //    throw NotSupported(
    //        "Not supported yet. Either ODENetworkModel or NetworkModel must be given.");
    //}
    model_ = model;
}

//void ODEWorld_New::bind_to(boost::shared_ptr<NetworkModel> model)
//{
//    if (boost::shared_ptr<NetworkModel> bound_model = model_.lock())
//    {
//        if (bound_model.get() != model.get())
//        {
//            std::cerr << "Warning: ODENetworkModel is already bound to ODEWorld."
//                << std::endl;
//        }
//    }
//    else if (generated_)
//    {
//        std::cerr << "Warning: NetworkModel is already bound to ODEWorld."
//            << std::endl;
//    }
//
//    this->model_ = model;
//    generated_.reset();
//}


void ODEWorld_New::save(const std::string& filename) const
{
    throw NotImplemented("ODEWorld_new::save not implemented");
}

void ODEWorld_New::load(const std::string& filename)
{
    throw NotImplemented("ODEWorld_new::load not implemented");
}


ODESimulator_New::deriv_func 
ODESimulator_New::generate_system()  const
{

    const std::vector<Species> species(world_->list_species());
    const reaction_rule_container_type& reaction_rules(model_->reaction_rules());

    // Step 1.  Assign index for each species.
    //   These indices are the order of species in dxdt container.
    typedef utils::get_mapper_mf<
        Species, state_type::size_type>::type species_map_type;
    species_map_type index_map;
    state_type::size_type i(0);   
    for(std::vector<Species>::const_iterator it(species.begin()); it != species.end(); it++) 
    {
        index_map[*it] = i; i++;
    }
    // Step 2. Create mapped_reaction_rule_container.
    mapped_reaction_container_type mapped_reactions;
    mapped_reactions.reserve(reaction_rules.size());
    for(reaction_rule_container_type::const_iterator it(reaction_rules.begin());
            it != reaction_rules.end(); it++) 
    {
        mapped_reaction_type mapped_rr;
        const ReactionRule::reactant_container_type &reactants(it->reactants());
        const ReactionRule::product_container_type  &products(it->products());
        for(ReactionRule::reactant_container_type::const_iterator r_it = reactants.begin(); 
                r_it != reactants.end(); r_it++) 
        {   
            mapped_rr.reactants.push_back(index_map[*r_it]);    
        }
        for(ReactionRule::product_container_type::const_iterator p_it = products.begin(); 
                p_it != products.end(); p_it++) 
        {
            mapped_rr.products.push_back(index_map[*p_it]);
        }
        mapped_rr.k = it->k();
        // Set default value for coefficients.
        //  These values may be over-written when reaction_rule_descriptor exist.
        mapped_rr.reactant_coefficients.resize( reactants.size() );
        mapped_rr.product_coefficients.resize( products.size() );
        std::fill(mapped_rr.reactant_coefficients.begin(), mapped_rr.reactant_coefficients.end(), 1);
        std::fill(mapped_rr.product_coefficients.begin(), mapped_rr.product_coefficients.end(), 1);
        // Additilnal Information(Overwrite default values defined above)
        if (it->has_descriptor()) {
            boost::shared_ptr<ReactionRuleDescriptor> rr_desc = it->get_descriptor();
            // Coefficients on both side
            if (0 < rr_desc->reactant_coefficients().size()) {
                mapped_rr.reactant_coefficients.clear();
                const ReactionRuleDescriptor::reaction_coefficient_list_type &reactant_coefficients(
                        rr_desc->reactant_coefficients());
                std::copy(reactant_coefficients.begin(), reactant_coefficients.end(), 
                        std::back_inserter(mapped_rr.reactant_coefficients) );
            }
            if (0 < rr_desc->product_coefficients().size()) {
                mapped_rr.product_coefficients.clear();
                const ReactionRuleDescriptor::reaction_coefficient_list_type &product_coefficients(
                        rr_desc->product_coefficients());
                std::copy(product_coefficients.begin(), product_coefficients.end(), 
                        std::back_inserter(mapped_rr.product_coefficients) );
            }
            // Reaction Propensity such as ratelaw
        } else {
            ;
        }
        mapped_reactions.push_back(mapped_rr);
    }
    // Step 3. Create Derivative and Jacobian object that bind mapped_reaction_container.
    return deriv_func(mapped_reactions, world_->volume() );
}

bool ODESimulator_New::step(const Real &upto)
{
    if (upto <= t()) {
        return false;
    }
    const Real dt(std::min(dt_,upto - t()));
    const Real ntime(std::min(upto, t() + dt_));
    const std::vector<Species> species(world_->list_species());
    state_type x(species.size());
    state_type::size_type i(0);
    // 1. Setup the array to create odeint object.
    for (NetworkModel::species_container_type::const_iterator it(species.begin());
            it != species.end(); it++) 
    {
        x[i] = static_cast<double>(world_->get_value_exact(*it));
        i++;
    }
    //std::pair<deriv_func, jacobi_func> system(generate_system());
    deriv_func  system(generate_system());
    StateAndTimeBackInserter::state_container_type x_vec;
    StateAndTimeBackInserter::time_container_type times;
    size_t steps;
    // 2. Instantiate the odeint object and run.
    switch (this->solver_type_) {
        case ecell4::ode::RUNGE_KUTTA_CASH_KARP54:
            {
                // This solver doesn't need the jacobian
                typedef odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;
                steps = (
                    odeint::integrate_adaptive(odeint::make_controlled<error_stepper_type>(abs_tol_, rel_tol_), 
                        system, x, t(), ntime, dt, StateAndTimeBackInserter(x_vec, times)));
                        //system.first, x, t(), ntime, dt, StateAndTimeBackInserter(x_vec, times)));
            }
            break;
        case ecell4::ode::ROSENBROCK4_CONTROLLER:
            {
                throw NotSupported("Sorry! Ronsenbrock4 is not implemented yet");
            }
            break;
        case ecell4::ode::EULER:
            {
                typedef odeint::euler<state_type> stepper_type;
                steps = (
                    odeint::integrate_const(
                        //stepper_type(), system.first, x, t(), ntime, dt,
                        stepper_type(), system, x, t(), ntime, dt,
                        StateAndTimeBackInserter(x_vec, times)));
            }
            break;
        default:
            throw IllegalState("Solver is not specified\n");
    }

    // 3. Write back the result of numerical-integration into the world class.
    {
        state_type::size_type i(0);
        for(Model::species_container_type::const_iterator it(species.begin());
                it != species.end(); it++)
        {
            world_->set_value(*it, static_cast<Real>(x_vec[steps](i)));
            i++;
        }
    }
    set_t(ntime);
    num_steps_++;
    return (ntime < upto);
}

} // ode
} // ecell4
