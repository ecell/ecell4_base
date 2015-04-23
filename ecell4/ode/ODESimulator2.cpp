#include "ODESimulator2.hpp"

#include <boost/numeric/odeint.hpp>
namespace odeint = boost::numeric::odeint;

namespace ecell4
{

namespace ode
{
std::pair<ODESimulator2::deriv_func, ODESimulator2::jacobi_func>
ODESimulator2::generate_system() const
{
    const std::vector<Species> species(world_->list_species());
    const ODENetworkModel::ode_reaction_rule_container_type& ode_reaction_rules(model_->ode_reaction_rules());
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
    reactions.reserve(ode_reaction_rules.size());
    for(ODENetworkModel::ode_reaction_rule_container_type::const_iterator
        i(ode_reaction_rules.begin()); i != ode_reaction_rules.end(); i++)
    {
        
    }
    return std::make_pair(
            deriv_func(reactions, world_->volume()),
            jacobi_func(reactions, world_->volume()));
}

bool ODESimulator2::step(const Real &upto)
{
    if (upto <= t())
    {
        return false;
    }
    const Real dt(upto - t());

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

    typedef odeint::rosenbrock4<state_type::value_type> error_stepper_type;
    typedef odeint::rosenbrock4_controller<error_stepper_type>
        controlled_stepper_type;

    const double abs_err(1e-10), rel_err(1e-6), a_x(1.0), a_dxdt(1.0);
    controlled_stepper_type controlled_stepper(abs_err, rel_err);
    const size_t steps(
        odeint::integrate_adaptive(
            controlled_stepper, system, x, t(), upto, dt,
            StateAndTimeBackInserter(x_vec, times)));
    {
        state_type::size_type i(0);
        for(ODENetworkModel::species_container_type::const_iterator
            it(species.begin()); it != species.end(); it++)
        {
            world_->set_value(*it, static_cast<Real>(x_vec[steps](i)));
            i++;
        }
    }
    set_t(upto);
    num_steps_++;
    return false;
}

} // ode
} // ecell4
