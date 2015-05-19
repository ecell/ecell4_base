#ifndef __ECELL4_ODE_ODE_SIMULATOR2_HPP
#define __ECELL4_ODE_ODE_SIMULATOR2_HPP

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

#include "ODEWorld.hpp"
#include "ODEReactionRule.hpp"
#include "ODENetworkModel.hpp"

namespace ecell4
{
namespace ode
{

class ODESimulator2
{
public:
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
        const Real volume_;
        const Real vinv_;
        const reaction_container_type reactions_;
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
                // prepare the matrix object that will be filled by ODERatelaw.
                matrix_type::size_type row_length = reactants_size + products_size;
                matrix_type::size_type col_length = row_length;
                matrix_type mat(row_length, col_length);
                
                // Call the ODERatelaw object
                if (i->ratelaw.expired() || i->ratelaw.lock()->is_available() == false)
                {
                    boost::scoped_ptr<ODERatelaw> temporary_ratelaw_obj(new ODERatelawMassAction(i->k));
                    temporary_ratelaw_obj->jacobi_func(mat, reactants_states, products_states, volume_, t, *(i->raw) );
                }
                else
                {
                    boost::shared_ptr<ODERatelaw> ratelaw = i->ratelaw.lock();
                    ratelaw->jacobi_func(mat, reactants_states, products_states, volume_, t, *(i->raw) );
                }

                // Merge matrix into whole system's jacobian.
                for(matrix_type::size_type row(0); row < row_length; row++)
                {
                    matrix_type::size_type j_row(row < reactants_size ? i->reactants[row] : i->products[row - reactants_size]);
                    for(matrix_type::size_type col(0); col < col_length; col++)
                    {
                        matrix_type::size_type j_col(col < reactants_size ? i->reactants[col] : i->products[col - reactants_size]);
                        jacobi(j_row, j_col) += mat(row, col);
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
    ODESimulator2(
        const boost::shared_ptr<ODENetworkModel> &model,
        const boost::shared_ptr<ODEWorld> &world)
        : model_(model), world_(world), dt_(inf), num_steps_(0)
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
    }
    bool step(const Real &upto);

    Real next_time() const
    {
        return this->t() + this->dt();
    }
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
    Integer num_steps() const
    {
        return this->num_steps_;
    }
protected:
    std::pair<deriv_func, jacobi_func> generate_system() const;
protected:
    Real dt_;
    boost::shared_ptr<ODENetworkModel> model_;
    boost::shared_ptr<ODEWorld> world_;

    Integer num_steps_;
};

} // ode

} // ecell4

#endif  // __ECELL4_ODE_ODE_SIMULATOR2_HPP
