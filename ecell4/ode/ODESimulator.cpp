#include "ODESimulator.hpp"

#include <boost/numeric/odeint.hpp>
namespace odeint = boost::numeric::odeint;

namespace ecell4
{

namespace ode
{

ODESystem ODESimulator::generate_system() const
{
    const std::vector<Species> species(world_->list_species());
    const Model::reaction_rule_container_type& reaction_rules(model_->reaction_rules());

    typedef utils::get_mapper_mf<
        Species, ODESystem::state_type::size_type>::type species_map_type;

    species_map_type index_map;
    ODESystem::state_type::size_type i(0);
    for (std::vector<Species>::const_iterator
        it(species.begin()); it != species.end(); ++it)
    {
        index_map[*it] = i;
        ++i;
    }

    std::vector<ODESystem::reaction_type> reactions;
    reactions.reserve(reaction_rules.size());
    for (Model::reaction_rule_container_type::const_iterator
        i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
    {
        const ReactionRule::reactant_container_type&
            reactants((*i).reactants());
        const ReactionRule::product_container_type&
            products((*i).products());

        ODESystem::reaction_type r;
        r.k = (*i).k();
        r.reactants.reserve(reactants.size());
        r.products.reserve(products.size());

        for (ReactionRule::reactant_container_type::const_iterator
                 j(reactants.begin()); j != reactants.end(); ++j)
        {
            r.reactants.push_back(index_map[*j]);
        }

        for (ReactionRule::product_container_type::const_iterator
                 j(products.begin()); j != products.end(); ++j)
        {
            r.products.push_back(index_map[*j]);
        }

        reactions.push_back(r);
    }

    return ODESystem(reactions, world_->volume());
}

bool ODESimulator::step(const Real& upto)
{
    if (upto <= t())
    {
        return false;
    }

    // initialize();

    const std::vector<Species> species(world_->list_species());

    ODESystem::state_type x(species.size());
    ODESystem::state_type::size_type i(0);
    for (Model::species_container_type::const_iterator
        it(species.begin()); it != species.end(); ++it)
    {
        x[i] = static_cast<double>(world_->get_value(*it));
        ++i;
    }

    typedef odeint::runge_kutta_cash_karp54<ODESystem::state_type>
        error_stepper_type;
    typedef odeint::controlled_runge_kutta<error_stepper_type>
        controlled_stepper_type;

    ODESystem func_obj(generate_system());
    StateAndTimeBackInserter::state_container_type x_vec;
    StateAndTimeBackInserter::time_container_type times;

    const Real dt(upto - t());

    // size_t steps(odeint::integrate(
    //                  func_obj, x, t(), upto, dt,
    //                  StateAndTimeBackInserter(x_vec, times)));

    const double abs_err(1e-10), rel_err(1e-6), a_x(1.0), a_dxdt(1.0);
    controlled_stepper_type controlled_stepper(
        odeint::default_error_checker<
            double, odeint::range_algebra, odeint::default_operations>(
                abs_err, rel_err, a_x, a_dxdt));
    const size_t steps(
        odeint::integrate_adaptive(
            controlled_stepper, func_obj, x, t(), upto, dt,
            StateAndTimeBackInserter(x_vec, times)));

    {
        ODESystem::state_type::size_type i(0);
        for (NetworkModel::species_container_type::const_iterator
                 it(species.begin()); it != species.end(); ++it)
        {
            world_->set_value(*it, static_cast<Real>(x_vec[steps][i]));
            ++i;
        }
    }

    set_t(upto);
    // set_t(times[steps]);
    num_steps_++;
    return false;
}

} // ode

} // ecell4
