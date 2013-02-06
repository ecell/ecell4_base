#include "ODESimulator.hpp"

#include <boost/numeric/odeint.hpp>
namespace odeint = boost::numeric::odeint;

namespace ecell4
{

namespace ode
{

bool ODESimulator::step(Real const& upto)
{
    NetworkModel::species_container_type const& species(model_->species());
    ODESystem::state_type x(species.size());

    {
        ODESystem::state_type::size_type i(0);
        for (NetworkModel::species_container_type::const_iterator
                 it(species.begin()); it != species.end(); ++it)
        {
            x[i] = static_cast<double>(world_->num_molecules(*it));
            ++i;
        }
    }

    typedef odeint::runge_kutta_cash_karp54<ODESystem::state_type>
        error_stepper_type;
    typedef odeint::controlled_runge_kutta<error_stepper_type>
        controlled_stepper_type;

    ODESystem func_obj(model_, world_->volume());
    StateAndTimeBackInserter::state_container_type x_vec;
    StateAndTimeBackInserter::time_container_type times;

    // size_t steps(odeint::integrate(
    //                  func_obj, x, t_, upto, upto - t_,
    //                  StateAndTimeBackInserter(x_vec, times)));

    double abs_err(1e-10), rel_err(1e-6), a_x(1.0), a_dxdt(1.0);
    controlled_stepper_type controlled_stepper(
        odeint::default_error_checker<double>(abs_err, rel_err, a_x, a_dxdt));
    size_t steps(odeint::integrate_adaptive(
                     controlled_stepper, func_obj, x, t_, upto, upto - t_,
                     StateAndTimeBackInserter(x_vec, times)));

    {
        ODESystem::state_type::size_type i(0);
        for (NetworkModel::species_container_type::const_iterator
                 it(species.begin()); it != species.end(); ++it)
        {
            world_->set_num_molecules(*it, static_cast<Real>(x_vec[steps][i]));
            ++i;
        }
    }

    t_ = upto;
    // t_ = times[steps];
    ++num_steps_;
    return false;
}

} // ode

} // ecell4
