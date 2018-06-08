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
    for(Model::reaction_rule_container_type::const_iterator
        i(reaction_rules.begin()); i != reaction_rules.end(); i++)
    {
        const ReactionRule& rr = (*i);

        const ReactionRule::reactant_container_type reactants(rr.reactants());
        const ReactionRule::product_container_type products(rr.products());

        reaction_type r;

        r.k = rr.k();

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

            const ReactionRuleDescriptor::coefficient_container_type& reactants_coeff = rrd->reactant_coefficients();
            std::copy(reactants_coeff.begin(), reactants_coeff.end(), std::back_inserter(r.reactant_coefficients));

            const ReactionRuleDescriptor::coefficient_container_type& products_coeff = rrd->product_coefficients();
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
