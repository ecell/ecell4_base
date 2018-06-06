#include "ODESimulator_New.hpp"

#include <boost/numeric/odeint.hpp>
#include <algorithm>

namespace odeint = boost::numeric::odeint;

namespace ecell4
{

namespace ode
{

std::pair<ODESimulator_New::deriv_func, ODESimulator_New::jacobi_func>
ODESimulator_New::generate_system() const
{
    const std::vector<Species> species(world_->list_species());
    const ODENetworkModel::ode_reaction_rule_container_type& ode_reaction_rules = ode_reaction_rules_;
    // const ODENetworkModel::ode_reaction_rule_container_type& ode_reaction_rules(model_->ode_reaction_rules());
    const Model::reaction_rule_container_type& reaction_rules = model_->reaction_rules();
    typedef utils::get_mapper_mf<
        Species, state_type::size_type>::type species_map_type;

    species_map_type index_map;
    state_type::size_type i(0);
    for(std::vector<Species>::const_iterator it(species.begin());
        it != species.end(); it++)
    {
        index_map[*it] = i;
        i++;
    }

    std::vector<reaction_type> reactions;
    reactions.reserve(reaction_rules.size());
    // reactions.reserve(ode_reaction_rules.size());
    assert(reaction_rules.size() == ode_reaction_rules.size());
    for (size_t i = 0; i < reaction_rules.size(); ++i)
    // for(ODENetworkModel::ode_reaction_rule_container_type::const_iterator
    //     i(ode_reaction_rules.begin()); i != ode_reaction_rules.end(); i++)
    {
        const ReactionRule& rr = reaction_rules[i];
        const ODEReactionRule& orr = ode_reaction_rules[i];

        const ReactionRule::reactant_container_type reactants(rr.reactants());
        const ReactionRule::product_container_type products(rr.products());

        reaction_type r;

        // r.raw = &(orr);
        r.k = rr.k();
        // if (orr.has_ratelaw())
        // {
        //     r.ratelaw = orr.get_ratelaw();
        // }

        r.reactants.reserve(reactants.size());
        for(ReactionRule::reactant_container_type::const_iterator j(reactants.begin());
            j != reactants.end(); j++)
        {
            r.reactants.push_back(index_map[*j]);
        }

        r.products.reserve(products.size());
        for(ReactionRule::product_container_type::const_iterator j(products.begin());
            j != products.end(); j++)
        {
            r.products.push_back(index_map[*j]);
        }

        if (rr.has_descriptor() && rr.get_descriptor()->has_coefficients())
        {
            const boost::shared_ptr<ReactionRuleDescriptor>& rrd(rr.get_descriptor());

            r.ratelaw = rr.get_descriptor();

            const ReactionRuleDescriptor::reaction_coefficient_list_type& reactants_coeff = rrd->reactant_coefficients();
            std::copy(reactants_coeff.begin(), reactants_coeff.end(), std::back_inserter(r.reactant_coefficients));

            const ReactionRuleDescriptor::reaction_coefficient_list_type& products_coeff = rrd->product_coefficients();
            std::copy(products_coeff.begin(), products_coeff.end(), std::back_inserter(r.product_coefficients));
        }
        else
        {
            r.reactant_coefficients.resize(reactants.size(), 1.0);
            r.product_coefficients.resize(products.size(), 1.0);
        }

        reactions.push_back(r);
    }
    return std::make_pair(
            deriv_func(reactions, world_->volume()),
            jacobi_func(reactions, world_->volume(), abs_tol_, rel_tol_));
}

bool ODESimulator_New::step(const Real &upto)
{
    if (upto <= t())
    {
        return false;
    }
    const Real dt(std::min(dt_, upto - t()));

    const Real ntime(std::min(upto, t() + dt_));

    //initialize();
    const std::vector<Species> species(world_->list_species());

    state_type x(species.size());
    state_type::size_type i(0);
    for(ODENetworkModel::species_container_type::const_iterator it(species.begin());
        it != species.end(); it++)
    {
        x[i] = static_cast<double>(world_->get_value_exact(*it));
        i++;
    }
    std::pair<deriv_func, jacobi_func> system(generate_system());
    StateAndTimeBackInserter::state_container_type x_vec;
    StateAndTimeBackInserter::time_container_type times;

    size_t steps;
    switch (this->solver_type_) {
        case ecell4::ode::RUNGE_KUTTA_CASH_KARP54:
            {
                /* This solver doesn't need the jacobian */
                typedef odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;
                steps = (
                    odeint::integrate_adaptive(
                        odeint::make_controlled<error_stepper_type>(abs_tol_, rel_tol_),
                        system.first, x, t(), ntime, dt,
                        StateAndTimeBackInserter(x_vec, times)));
            }
            break;
        case ecell4::ode::ROSENBROCK4_CONTROLLER:
            {
                typedef odeint::rosenbrock4<state_type::value_type> error_stepper_type;
                steps = (
                    odeint::integrate_adaptive(
                        odeint::make_controlled<error_stepper_type>(abs_tol_, rel_tol_),
                        system, x, t(), ntime, dt,
                        StateAndTimeBackInserter(x_vec, times)));
            }
            break;
        case ecell4::ode::EULER:
            {
                typedef odeint::euler<state_type> stepper_type;
                steps = (
                    odeint::integrate_const(
                        stepper_type(), system.first, x, t(), ntime, dt,
                        StateAndTimeBackInserter(x_vec, times)));
            }
            break;
        default:
            throw IllegalState("Solver is not specified\n");
    };

    // const double a_x(1.0), a_dxdt(1.0);
    // const size_t steps(
    //     odeint::integrate_adaptive(
    //         controlled_stepper, system, x, t(), upto, dt_,
    //         StateAndTimeBackInserter(x_vec, times)));
    {
        state_type::size_type i(0);
        for(ODENetworkModel::species_container_type::const_iterator
            it(species.begin()); it != species.end(); it++)
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

//XXX: #include "ODESimulator_New.hpp"
//XXX: 
//XXX: #include <boost/numeric/odeint.hpp>
//XXX: #include <algorithm>
//XXX: #include <functional>
//XXX: 
//XXX: namespace odeint = boost::numeric::odeint;
//XXX: 
//XXX: namespace ecell4
//XXX: {
//XXX: 
//XXX: namespace ode
//XXX: {
//XXX: 
//XXX: 
//XXX: void ODEWorld_New::bind_to(boost::shared_ptr<Model> model)
//XXX: {
//XXX:     //if (generated_)
//XXX:     //{
//XXX:     //    std::cerr << "Warning: NetworkModel is already bound to ODEWorld."
//XXX:     //        << std::endl;
//XXX:     //}
//XXX:     //else if (model_.expired())
//XXX:     //{
//XXX:     //    std::cerr << "Warning: ODENetworkModel is already bound to ODEWorld."
//XXX:     //        << std::endl;
//XXX:     //}
//XXX:     //try
//XXX:     //{
//XXX:     //    boost::shared_ptr<NetworkModel> tmp(new NetworkModel(model));
//XXX:     //    generated_.swap(tmp);
//XXX:     //    model_.reset();
//XXX:     //}
//XXX:     //catch (NotSupported e)
//XXX:     //{
//XXX:     //    throw NotSupported(
//XXX:     //        "Not supported yet. Either ODENetworkModel or NetworkModel must be given.");
//XXX:     //}
//XXX:     model_ = model;
//XXX: }
//XXX: 
//XXX: //void ODEWorld_New::bind_to(boost::shared_ptr<NetworkModel> model)
//XXX: //{
//XXX: //    if (boost::shared_ptr<NetworkModel> bound_model = model_.lock())
//XXX: //    {
//XXX: //        if (bound_model.get() != model.get())
//XXX: //        {
//XXX: //            std::cerr << "Warning: ODENetworkModel is already bound to ODEWorld."
//XXX: //                << std::endl;
//XXX: //        }
//XXX: //    }
//XXX: //    else if (generated_)
//XXX: //    {
//XXX: //        std::cerr << "Warning: NetworkModel is already bound to ODEWorld."
//XXX: //            << std::endl;
//XXX: //    }
//XXX: //
//XXX: //    this->model_ = model;
//XXX: //    generated_.reset();
//XXX: //}
//XXX: 
//XXX: 
//XXX: void ODEWorld_New::save(const std::string& filename) const
//XXX: {
//XXX:     throw NotImplemented("ODEWorld_new::save not implemented");
//XXX: }
//XXX: 
//XXX: void ODEWorld_New::load(const std::string& filename)
//XXX: {
//XXX:     throw NotImplemented("ODEWorld_new::load not implemented");
//XXX: }
//XXX: 
//XXX: 
//XXX: ODESimulator_New::deriv_func 
//XXX: ODESimulator_New::generate_system()  const
//XXX: {
//XXX: 
//XXX:     const std::vector<Species> species(world_->list_species());
//XXX:     const reaction_rule_container_type& reaction_rules(model_->reaction_rules());
//XXX: 
//XXX:     // Step 1.  Assign index for each species.
//XXX:     //   These indices are the order of species in dxdt container.
//XXX:     typedef utils::get_mapper_mf<
//XXX:         Species, state_type::size_type>::type species_map_type;
//XXX:     species_map_type index_map;
//XXX:     state_type::size_type i(0);   
//XXX:     for(std::vector<Species>::const_iterator it(species.begin()); it != species.end(); it++) 
//XXX:     {
//XXX:         index_map[*it] = i; i++;
//XXX:     }
//XXX:     // Step 2. Create mapped_reaction_rule_container.
//XXX:     mapped_reaction_container_type mapped_reactions;
//XXX:     mapped_reactions.reserve(reaction_rules.size());
//XXX:     for(reaction_rule_container_type::const_iterator it(reaction_rules.begin());
//XXX:             it != reaction_rules.end(); it++) 
//XXX:     {
//XXX:         mapped_reaction_type mapped_rr;
//XXX:         const ReactionRule::reactant_container_type &reactants(it->reactants());
//XXX:         const ReactionRule::product_container_type  &products(it->products());
//XXX:         for(ReactionRule::reactant_container_type::const_iterator r_it = reactants.begin(); 
//XXX:                 r_it != reactants.end(); r_it++) 
//XXX:         {   
//XXX:             mapped_rr.reactants.push_back(index_map[*r_it]);    
//XXX:         }
//XXX:         for(ReactionRule::product_container_type::const_iterator p_it = products.begin(); 
//XXX:                 p_it != products.end(); p_it++) 
//XXX:         {
//XXX:             mapped_rr.products.push_back(index_map[*p_it]);
//XXX:         }
//XXX:         mapped_rr.k = it->k();
//XXX:         // Set default value for coefficients.
//XXX:         //  These values may be over-written when reaction_rule_descriptor exist.
//XXX:         mapped_rr.reactant_coefficients.resize( reactants.size() );
//XXX:         mapped_rr.product_coefficients.resize( products.size() );
//XXX:         std::fill(mapped_rr.reactant_coefficients.begin(), mapped_rr.reactant_coefficients.end(), 1);
//XXX:         std::fill(mapped_rr.product_coefficients.begin(), mapped_rr.product_coefficients.end(), 1);
//XXX:         // Additilnal Information(Overwrite default values defined above)
//XXX:         if (it->has_descriptor()) {
//XXX:             boost::shared_ptr<ReactionRuleDescriptor> rr_desc = it->get_descriptor();
//XXX:             // Coefficients on both side
//XXX:             if (0 < rr_desc->reactant_coefficients().size()) {
//XXX:                 mapped_rr.reactant_coefficients.clear();
//XXX:                 const ReactionRuleDescriptor::reaction_coefficient_list_type &reactant_coefficients(
//XXX:                         rr_desc->reactant_coefficients());
//XXX:                 std::copy(reactant_coefficients.begin(), reactant_coefficients.end(), 
//XXX:                         std::back_inserter(mapped_rr.reactant_coefficients) );
//XXX:             }
//XXX:             if (0 < rr_desc->product_coefficients().size()) {
//XXX:                 mapped_rr.product_coefficients.clear();
//XXX:                 const ReactionRuleDescriptor::reaction_coefficient_list_type &product_coefficients(
//XXX:                         rr_desc->product_coefficients());
//XXX:                 std::copy(product_coefficients.begin(), product_coefficients.end(), 
//XXX:                         std::back_inserter(mapped_rr.product_coefficients) );
//XXX:             }
//XXX:             // Reaction Propensity such as ratelaw
//XXX:         } else {
//XXX:             ;
//XXX:         }
//XXX:         mapped_reactions.push_back(mapped_rr);
//XXX:     }
//XXX:     // Step 3. Create Derivative and Jacobian object that bind mapped_reaction_container.
//XXX:     return deriv_func(mapped_reactions, world_->volume() );
//XXX: }
//XXX: 
//XXX: bool ODESimulator_New::step(const Real &upto)
//XXX: {
//XXX:     if (upto <= t()) {
//XXX:         return false;
//XXX:     }
//XXX:     const Real dt(std::min(dt_,upto - t()));
//XXX:     const Real ntime(std::min(upto, t() + dt_));
//XXX:     const std::vector<Species> species(world_->list_species());
//XXX:     state_type x(species.size());
//XXX:     state_type::size_type i(0);
//XXX:     // 1. Setup the array to create odeint object.
//XXX:     for (NetworkModel::species_container_type::const_iterator it(species.begin());
//XXX:             it != species.end(); it++) 
//XXX:     {
//XXX:         x[i] = static_cast<double>(world_->get_value_exact(*it));
//XXX:         i++;
//XXX:     }
//XXX:     //std::pair<deriv_func, jacobi_func> system(generate_system());
//XXX:     deriv_func  system(generate_system());
//XXX:     StateAndTimeBackInserter::state_container_type x_vec;
//XXX:     StateAndTimeBackInserter::time_container_type times;
//XXX:     size_t steps;
//XXX:     // 2. Instantiate the odeint object and run.
//XXX:     switch (this->solver_type_) {
//XXX:         case ecell4::ode::RUNGE_KUTTA_CASH_KARP54:
//XXX:             {
//XXX:                 // This solver doesn't need the jacobian
//XXX:                 typedef odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;
//XXX:                 steps = (
//XXX:                     odeint::integrate_adaptive(odeint::make_controlled<error_stepper_type>(abs_tol_, rel_tol_), 
//XXX:                         system, x, t(), ntime, dt, StateAndTimeBackInserter(x_vec, times)));
//XXX:                         //system.first, x, t(), ntime, dt, StateAndTimeBackInserter(x_vec, times)));
//XXX:             }
//XXX:             break;
//XXX:         case ecell4::ode::ROSENBROCK4_CONTROLLER:
//XXX:             {
//XXX:                 throw NotSupported("Sorry! Ronsenbrock4 is not implemented yet");
//XXX:             }
//XXX:             break;
//XXX:         case ecell4::ode::EULER:
//XXX:             {
//XXX:                 typedef odeint::euler<state_type> stepper_type;
//XXX:                 steps = (
//XXX:                     odeint::integrate_const(
//XXX:                         //stepper_type(), system.first, x, t(), ntime, dt,
//XXX:                         stepper_type(), system, x, t(), ntime, dt,
//XXX:                         StateAndTimeBackInserter(x_vec, times)));
//XXX:             }
//XXX:             break;
//XXX:         default:
//XXX:             throw IllegalState("Solver is not specified\n");
//XXX:     }
//XXX: 
//XXX:     // 3. Write back the result of numerical-integration into the world class.
//XXX:     {
//XXX:         state_type::size_type i(0);
//XXX:         for(Model::species_container_type::const_iterator it(species.begin());
//XXX:                 it != species.end(); it++)
//XXX:         {
//XXX:             world_->set_value(*it, static_cast<Real>(x_vec[steps](i)));
//XXX:             i++;
//XXX:         }
//XXX:     }
//XXX:     set_t(ntime);
//XXX:     num_steps_++;
//XXX:     return (ntime < upto);
//XXX: }
//XXX: 
//XXX: } // ode
//XXX: } // ecell4
