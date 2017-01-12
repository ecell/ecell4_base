#include <ecell4/core/exceptions.hpp>

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "ODERatelaw.hpp"
#include "ODEReactionRule.hpp"

namespace ecell4
{
namespace ode 
{
Real ODERatelawCppCallback::deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, 
        Real const volume, Real const t,
        ODEReactionRule const &rr)
{
    if (!is_available())
    {
        throw IllegalState("Callback Function has not been registerd");
    }
    return this->func_(reactants_state_array, products_state_array, volume, t, rr);
}

Real ODERatelawMassAction::deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, 
        Real const volume, Real const t,
        ODEReactionRule const &rr)
{
    ODEReactionRule::coefficient_container_type const reactants_coefficients(rr.reactants_coefficients());
    Real flux(this->k_ * volume);
    int i = 0;
    for(state_container_type::const_iterator it(reactants_state_array.begin());
        it != reactants_state_array.end(); it++, i++)
    {
        flux *= std::pow( (*it) / volume, reactants_coefficients[i]);
    }
    return flux;
}

Real ODERatelawCythonCallback::deriv_func(
    state_container_type const &reactants_state_array,
    state_container_type const &products_state_array, 
    Real const volume, Real const t,
    ODEReactionRule const &rr)
{
    if (!is_available())
    {
        throw IllegalState("Callback Function has not been registerd");
    }
    ODEReactionRule rr_tempolrary(rr);
    return this->indirect_func_(
        this->python_func_, reactants_state_array, products_state_array, 
        volume, t, &rr_tempolrary);
}

boost::shared_ptr<ODERatelawMassAction> to_ODERatelawMassAction(boost::shared_ptr<ODERatelaw> p)
{
    return boost::dynamic_pointer_cast<ODERatelawMassAction>(p);
}

boost::shared_ptr<ODERatelawCythonCallback> to_ODERatelawCythonCallback(boost::shared_ptr<ODERatelaw> p)
{
    return boost::dynamic_pointer_cast<ODERatelawCythonCallback>(p);
}

} // ode

} // ecell4
