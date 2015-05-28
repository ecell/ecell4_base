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
    return this->func_(reactants_state_array, products_state_array, volume);
}

void ODERatelawCppCallback::jacobi_func(
        matrix_type &jacobian,
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array,
        Real const volume, Real const t, 
        ODEReactionRule const &rr)
{
    Real h(this->h_); //XXX: 1.0e-8. Should be fixed

    std::fill(jacobian.data().begin(), jacobian.data().end(), Real(0.0));
    Real flux(this->deriv_func(reactants_state_array, products_state_array, volume, t, rr));
    double num_reactants(reactants_state_array.size());
    double num_products(products_state_array.size());

    // Differentiates by Reactants.
    for (int i(0); i < num_reactants; i++)
    {
        //XXX: For now, we are using FORWARD difference method.
        state_container_type h_shift(reactants_state_array);
        h_shift[i] += h;
        double deriv = (
            this->deriv_func(h_shift, products_state_array, volume, t, rr) - flux) / h;
        for (matrix_type::size_type j(0); j < jacobian.size1(); j++)
        {
            if (j < num_reactants)
            {
                jacobian(j, i) -= deriv;
            }
            else
            {
                jacobian(j, i) += deriv;
            }
        }
    }

    // Differentiates by Products.
    for (int i(0); i < num_products; i++)
    {
        state_container_type h_shift(products_state_array);
        h_shift[i] += h;
        double deriv = (
            this->deriv_func(reactants_state_array, h_shift, volume, t, rr) - flux) / h;
        for (matrix_type::size_type j(0); j < jacobian.size1(); j++)
        {
            if (j < num_reactants)
            {
                jacobian(j, i + num_reactants) -= deriv;
            }
            else
            {
                jacobian(j, i + num_reactants) += deriv;
            }
        }
    }
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

void ODERatelawMassAction::jacobi_func(
        matrix_type &jacobian,
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, 
        Real const volume, Real const t, 
        ODEReactionRule const &rr)
{
    // The size of the argument 'state_array' must be resized
    // to the number of reactants.
    // The size of the argument 'jacobian' must be resized
    // to the number of (reactants + products)
    std::fill(jacobian.data().begin(), jacobian.data().end(), Real(0.0));
    Real flux(this->deriv_func(reactants_state_array, products_state_array, volume, t, rr));
    if (flux == Real(0.0))
    {
        return;
    }

    matrix_type::size_type num_reactants(
        static_cast<matrix_type::size_type>(reactants_state_array.size()));
    for(matrix_type::size_type i(0); i < num_reactants; i++)
    {
        Real partial(flux / reactants_state_array[i]);
        for(matrix_type::size_type j(0); j < jacobian.size1(); j++)
        {
            if (j < num_reactants)
            {
                jacobian(j, i) -= partial;
            }
            else
            {
                jacobian(j, i) += partial;
            }
        }
    }
}

} // ode

} // ecell4
