#ifndef ECELL4_ODE_ODE_SIMULATOR2_HPP
#define ECELL4_ODE_ODE_SIMULATOR2_HPP

#include <cstring>
#include <vector>
#include <numeric>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/get_mapper_mf.hpp>

#include <ecell4/core/SimulatorBase.hpp>

#include "ODEWorld.hpp"
#include "ODEReactionRule.hpp"
#include "ODENetworkModel.hpp"

namespace ecell4
{
namespace ode
{

enum ODESolverType {
    RUNGE_KUTTA_CASH_KARP54 = 0,
    ROSENBROCK4_CONTROLLER = 1,
    EULER = 2,
};

class ODESimulator
    : public SimulatorBase<ODENetworkModel, ODEWorld>
{
public:

    typedef SimulatorBase<ODENetworkModel, ODEWorld> base_type;

public:

    typedef boost::numeric::ublas::vector<double> state_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef std::vector<state_type::size_type> index_container_type;
    typedef std::vector<Real> coefficient_container_type;

    typedef ODEReactionRule reacton_container_type;

    struct reaction_type
    {
        index_container_type reactants;
        coefficient_container_type reactant_coefficients;
        index_container_type products;
        coefficient_container_type product_coefficients;
        Real k;
        boost::weak_ptr<ODERatelaw> ratelaw;
        const ODEReactionRule *raw;
    };
    typedef std::vector<reaction_type> reaction_container_type;

    class deriv_func
    {
    public:
        deriv_func(const reaction_container_type &reactions, const Real &volume)
            : reactions_(reactions), volume_(volume), vinv_(1.0 / volume)
        {
            ;
        }

        void operator()(const state_type &x, state_type &dxdt, const double &t)
        {
            std::fill(dxdt.begin(), dxdt.end(), 0.0);
            for(reaction_container_type::const_iterator i(reactions_.begin());
                i != reactions_.end(); i++)
            {
                ODERatelaw::state_container_type reactants_states(i->reactants.size());
                ODERatelaw::state_container_type products_states(i->products.size());
                ODERatelaw::state_container_type::size_type cnt(0);

                for(index_container_type::const_iterator j(i->reactants.begin());
                    j != i->reactants.end(); j++, cnt++)
                {
                    reactants_states[cnt] = x[*j];
                }
                cnt = 0;
                for(index_container_type::const_iterator j(i->products.begin());
                    j != i->products.end(); j++, cnt++)
                {
                    products_states[cnt] = x[*j];
                }
                double flux;
                // Calculation! XXX
                if (i->ratelaw.expired() || i->ratelaw.lock()->is_available() == false)
                {
                    boost::scoped_ptr<ODERatelaw> temporary_ratelaw_obj(new ODERatelawMassAction(i->k));
                    flux = temporary_ratelaw_obj->deriv_func(reactants_states, products_states, volume_, t, *(i->raw) );
                }
                else
                {
                    boost::shared_ptr<ODERatelaw> ratelaw = i->ratelaw.lock();
                    flux = ratelaw->deriv_func(reactants_states, products_states, volume_, t, *(i->raw) );
                }
                // Merge each reaction's flux into whole dxdt
                std::size_t nth = 0;
                for(index_container_type::const_iterator j(i->reactants.begin());
                    j != i->reactants.end(); j++)
                {
                    dxdt[*j] -= (flux * (double)i->reactant_coefficients[nth]);
                    nth++;
                }
                nth = 0;
                for(index_container_type::const_iterator j(i->products.begin()); 
                    j != i->products.end(); j++)
                {
                    dxdt[*j] += (flux * (double)i->product_coefficients[nth]);
                    nth++;
                }
            }
            return;
        }
    protected:
        const reaction_container_type reactions_;
        const Real volume_;
        const Real vinv_;
    };

    class jacobi_func
    {
    public:
        jacobi_func(const reaction_container_type &reactions, const Real& volume)
            : reactions_(reactions), volume_(volume), vinv_(1.0 / volume)
        {
            ;
        }
        void operator()(
                const state_type& x, matrix_type& jacobi, const double &t, state_type &dfdt) const
        {
            //fill 0 into jacobi and dfdt
            std::fill(dfdt.begin(), dfdt.end(), 0.0);
            std::fill(jacobi.data().begin(), jacobi.data().end(), 0.0);

            const Real h(1.0e-8);
            const Real ht(1.0e-10);

            // calculate jacobian for each reaction and merge it.
            for(reaction_container_type::const_iterator i(reactions_.begin()); 
                i != reactions_.end(); i++)
            {
                // Calculate one reactions's jabobian
                //  Prepare the state_array to pass ODERatelaw.
                index_container_type::size_type reactants_size(i->reactants.size());
                index_container_type::size_type products_size(i->products.size());
                ODERatelaw::state_container_type reactants_states(reactants_size);
                ODERatelaw::state_container_type products_states(products_size);
                ODERatelaw::state_container_type::size_type cnt(0);
                for(index_container_type::const_iterator j(i->reactants.begin());
                    j != i->reactants.end(); j++, cnt++)
                {
                    reactants_states[cnt] = x[*j];
                }
                cnt = 0;
                for(index_container_type::const_iterator j(i->products.begin());
                    j != i->products.end(); j++, cnt++)
                {
                    products_states[cnt] = x[*j];
                }
                // Call the ODERatelaw object
                if (i->ratelaw.expired() || i->ratelaw.lock()->is_available() == false)
                {
                    boost::scoped_ptr<ODERatelaw> temporary_ratelaw_obj(new ODERatelawMassAction(i->k));
                    Real flux_0 = temporary_ratelaw_obj->deriv_func(reactants_states, products_states, volume_, t, *(i->raw) );
                    // Differentiate by time
                    {
                        Real flux = temporary_ratelaw_obj->deriv_func(reactants_states, products_states, volume_, t + ht, *(i->raw) );
                        Real flux_deriv = (flux - flux_0) / h;
                        if (flux_deriv != 0.0)
                        {
                            for(std::size_t k(0); k < i->reactants.size(); k++)
                            {
                                matrix_type::size_type row = i->reactants[k];
                                Real coeff = i->reactant_coefficients[k];
                                dfdt[row] -= coeff * flux_deriv;
                            }
                            for(std::size_t k(0); k < i->products.size(); k++)
                            {
                                matrix_type::size_type row = i->products[k];
                                Real coeff = i->product_coefficients[k];
                                dfdt[row] += coeff * flux_deriv;
                            }
                        }
                    }
                    // Differentiate by each Reactants
                    for(std::size_t j(0); j < reactants_states.size(); j++)
                    {
                        ODERatelaw::state_container_type h_shift(reactants_states);
                        h_shift[j] += h;
                        Real flux = temporary_ratelaw_obj->deriv_func(h_shift, products_states, volume_, t, *(i->raw) );
                        Real flux_deriv = (flux - flux_0) / h;
                        matrix_type::size_type col = i->reactants[j];
                        for(std::size_t k(0); k < i->reactants.size(); k++)
                        {
                            matrix_type::size_type row = i->reactants[k];
                            Real coeff = i->reactant_coefficients[k];
                            jacobi(row, col) -= coeff * flux_deriv;
                        }
                        for(std::size_t k(0); k < i->products.size(); k++)
                        {
                            matrix_type::size_type row = i->products[k];
                            Real coeff = i->product_coefficients[k];
                            jacobi(row, col) += coeff * flux_deriv;
                        }
                    }
                    // Differentiate by Products
                    for(std::size_t j(0); j < products_states.size(); j++)
                    {
                        ODERatelaw::state_container_type h_shift(products_states);
                        h_shift[j] += h;
                        Real flux = temporary_ratelaw_obj->deriv_func(reactants_states, h_shift, volume_, t, *(i->raw));
                        Real flux_deriv = (flux - flux_0) / h;
                        matrix_type::size_type col = i->products[j];
                        for(std::size_t k(0); k < i->reactants.size(); k++)
                        {
                            matrix_type::size_type row = i->reactants[k];
                            Real coeff = i->reactant_coefficients[k];
                            jacobi(row, col) -= coeff * flux_deriv;
                        }
                        for(std::size_t k(0); k < i->products.size(); k++)
                        {
                            matrix_type::size_type row = i->products[k];
                            Real coeff = i->product_coefficients[k];
                            jacobi(row, col) += coeff * flux_deriv;
                        }
                    }
                }
                else
                {
                    boost::shared_ptr<ODERatelaw> ratelaw = i->ratelaw.lock();
                    Real flux_0 = ratelaw->deriv_func(reactants_states, products_states, volume_, t, *(i->raw) );
                    // Differentiate by time
                    {
                        Real flux = ratelaw->deriv_func(reactants_states, products_states, volume_, t + ht, *(i->raw) );
                        Real flux_deriv = (flux - flux_0) / h;
                        if (flux_deriv != 0.0)
                        {
                            for(std::size_t k(0); k < i->reactants.size(); k++)
                            {
                                matrix_type::size_type row = i->reactants[k];
                                Real coeff = i->reactant_coefficients[k];
                                dfdt[row] -= coeff * flux_deriv;
                            }
                            for(std::size_t k(0); k < i->products.size(); k++)
                            {
                                matrix_type::size_type row = i->products[k];
                                Real coeff = i->product_coefficients[k];
                                dfdt[row] += coeff * flux_deriv;
                            }
                        }
                    }
                    // Differentiate by each Reactants
                    for(std::size_t j(0); j < reactants_states.size(); j++)
                    {
                        ODERatelaw::state_container_type h_shift(reactants_states);
                        h_shift[j] += h;
                        Real flux = ratelaw->deriv_func(h_shift, products_states, volume_, t, *(i->raw) );
                        Real flux_deriv = (flux - flux_0) / h;
                        matrix_type::size_type col = i->reactants[j];
                        for(std::size_t k(0); k < i->reactants.size(); k++)
                        {
                            matrix_type::size_type row = i->reactants[k];
                            Real coeff = i->reactant_coefficients[k];
                            jacobi(row, col) -= coeff * flux_deriv;
                        }
                        for(std::size_t k(0); k < i->products.size(); k++)
                        {
                            matrix_type::size_type row = i->products[k];
                            Real coeff = i->product_coefficients[k];
                            jacobi(row, col) += coeff * flux_deriv;
                        }
                    }
                    // Differentiate by Products
                    for(std::size_t j(0); j < products_states.size(); j++)
                    {
                        ODERatelaw::state_container_type h_shift(products_states);
                        h_shift[j] += h;
                        Real flux = ratelaw->deriv_func(reactants_states, h_shift, volume_, t, *(i->raw));
                        Real flux_deriv = (flux - flux_0) / h;
                        matrix_type::size_type col = i->products[j];
                        for(std::size_t k(0); k < i->reactants.size(); k++)
                        {
                            matrix_type::size_type row = i->reactants[k];
                            Real coeff = i->reactant_coefficients[k];
                            jacobi(row, col) -= coeff * flux_deriv;
                        }
                        for(std::size_t k(0); k < i->products.size(); k++)
                        {
                            matrix_type::size_type row = i->products[k];
                            Real coeff = i->product_coefficients[k];
                            jacobi(row, col) += coeff * flux_deriv;
                        }
                    }
                }
            }
        }
    protected:
        const reaction_container_type reactions_;
        const Real volume_;
        const Real vinv_;
    };


    struct StateAndTimeBackInserter
    {
        typedef std::vector<state_type> state_container_type;
        typedef std::vector<double> time_container_type;

        state_container_type &m_states;
        time_container_type &m_times;
        StateAndTimeBackInserter(
            state_container_type &states, time_container_type &times)
            : m_states(states), m_times(times)
        {
            ;
        }
        void operator()(const state_type &x, double t)
        {
            m_states.push_back(x);
            m_times.push_back(t);
        }
    };
public:

    ODESimulator(
        const boost::shared_ptr<ODENetworkModel>& model,
        const boost::shared_ptr<ODEWorld>& world,
        const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
        : base_type(model, world), dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6),
          solver_type_(solver_type)
    {
        initialize();
    }

    ODESimulator(
        const boost::shared_ptr<ODEWorld>& world,
        const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
        : base_type(world), dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6),
          solver_type_(solver_type)
    {
        initialize();
    }

    ODESimulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<ODEWorld>& world,
        const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
        : base_type(boost::shared_ptr<ODENetworkModel>(new ODENetworkModel(model)), world),
          dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6), solver_type_(solver_type)
    {
        initialize();
    }

    void initialize()
    {
        const std::vector<Species> species(model_->list_species());
        for(std::vector<Species>::const_iterator it = species.begin();
                it != species.end(); it++)
        {
            if (!(world_->has_species(*it)))
            {
                world_->reserve_species(*it);
            }
        }
    }

    void step(void)
    {
        step(next_time());
        if ( this->model_->has_network_model() )
        {
            this->model_->update_model();
        }
    }
    bool step(const Real &upto);

    // Real next_time() const
    // {
    //     return this->t() + this->dt();
    // }
    // SimulatorTraits

    Real t(void) const
    {
        return world_->t();
    }
    void set_t(const Real &t)
    {
        world_->set_t(t);
    }
    Real dt(void) const
    {
        return this->dt_;
    }
    void set_dt(const Real &dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }
    // Integer num_steps() const
    // {
    //     return this->num_steps_;
    // }
    //

    Real absolute_tolerance() const
    {
        return abs_tol_;
    }

    void set_absolute_tolerance(const Real abs_tol)
    {
        if (abs_tol < 0)
        {
            throw std::invalid_argument("A tolerance must be positive or zero.");
        }
        abs_tol_ = abs_tol;
    }

    Real relative_tolerance() const
    {
        return rel_tol_;
    }

    void set_relative_tolerance(const Real rel_tol)
    {
        if (rel_tol < 0)
        {
            throw std::invalid_argument("A tolerance must be positive or zero.");
        }
        rel_tol_ = rel_tol;
    }

protected:
    std::pair<deriv_func, jacobi_func> generate_system() const;
protected:
    // boost::shared_ptr<ODENetworkModel> model_;
    // boost::shared_ptr<ODEWorld> world_;
    Real dt_;
    // Integer num_steps_;
    Real abs_tol_, rel_tol_;
    ODESolverType solver_type_;
};

} // ode

} // ecell4

#endif  // ECELL4_ODE_ODE_SIMULATOR2_HPP
