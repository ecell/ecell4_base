#ifndef ECELL4_ODE_ODE_SIMULATOR_NEW_HPP
#define ECELL4_ODE_ODE_SIMULATOR_NEW_HPP

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

#include "ODEWorld_New.hpp"
#include "ODEReactionRule.hpp"
#include "ODENetworkModel.hpp"

#include "ODESimulator.hpp"

namespace ecell4
{
namespace ode
{

// enum ODESolverType {
//     RUNGE_KUTTA_CASH_KARP54 = 0,
//     ROSENBROCK4_CONTROLLER = 1,
//     EULER = 2,
// };

class ODESimulator_New
    : public SimulatorBase<ODENetworkModel, ODEWorld_New>
{
public:

    typedef SimulatorBase<ODENetworkModel, ODEWorld_New> base_type;

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
        jacobi_func(
            const reaction_container_type &reactions, const Real& volume,
            const Real& abs_tol, const Real& rel_tol)
            : reactions_(reactions), volume_(volume), vinv_(1.0 / volume), abs_tol_(abs_tol), rel_tol_(rel_tol)
        {
            ;
        }
        void operator()(
                const state_type& x, matrix_type& jacobi, const double &t, state_type &dfdt) const
        {
            //fill 0 into jacobi and dfdt
            std::fill(dfdt.begin(), dfdt.end(), 0.0);
            std::fill(jacobi.data().begin(), jacobi.data().end(), 0.0);

            // const Real ETA(2.2204460492503131e-16);
            const Real SQRTETA(1.4901161193847656e-08);
            const Real r0(1.0);
            // Real fac(0.0);
            // for (std::size_t k(0); k < dfdt.size(); ++k)
            // {
            //     const Real ewtk(atol + rtol * x[k]);
            //     fac = std::max(fac, dfdt[k] * ewtk);
            // }
            // const Real r0(1000.0 * h * ETA * dfdt.size() * fac);  //XXX: h means the step interval
            // {
            //     const Real ewtj(atol + rtol * x[j]);
            //     const Real dyj(std::max(SQRTETA * abs(x[j]), r0 * ewtj));
            // }

            // const Real h(1.0e-8);
            // const Real ht(1.0e-10);
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
                        Real flux_deriv = (flux - flux_0) / ht;
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
                        const Real ewt = abs_tol_ + rel_tol_ * abs(reactants_states[j]);
                        const Real h = std::max(SQRTETA * abs(reactants_states[j]), r0 * ewt);
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
                        const Real ewt = abs_tol_ + rel_tol_ * abs(products_states[j]);
                        const Real h = std::max(SQRTETA * abs(products_states[j]), r0 * ewt);
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
                        Real flux_deriv = (flux - flux_0) / ht;
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
                        const Real ewt = abs_tol_ + rel_tol_ * abs(reactants_states[j]);
                        const Real h = std::max(SQRTETA * abs(reactants_states[j]), r0 * ewt);
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
                        const Real ewt = abs_tol_ + rel_tol_ * abs(products_states[j]);
                        const Real h = std::max(SQRTETA * abs(products_states[j]), r0 * ewt);
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
        const Real abs_tol_, rel_tol_;
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

    ODESimulator_New(
        const boost::shared_ptr<ODENetworkModel>& model,
        const boost::shared_ptr<ODEWorld_New>& world,
        const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
        : base_type(model, world), dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6),
          solver_type_(solver_type)
    {
        initialize();
    }

    ODESimulator_New(
        const boost::shared_ptr<ODEWorld_New>& world,
        const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
        : base_type(world), dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6),
          solver_type_(solver_type)
    {
        initialize();
    }

    ODESimulator_New(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<ODEWorld_New>& world,
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

    Real evaluate(const ReactionRule& rr) const
    {
        return evaluate(ODEReactionRule(rr));
    }

    Real evaluate(const ODEReactionRule& rr) const
    {
        if (!rr.has_ratelaw())
        {
            // std::cout << "No ratelaw was bound. Return zero." << std::endl;
            return 0.0;
        }

        const ODEReactionRule::reactant_container_type reactants = rr.reactants();
        const ODEReactionRule::product_container_type products = rr.products();

        ODERatelaw::state_container_type::size_type cnt(0);

        ODERatelaw::state_container_type r(reactants.size());
        for(ODEReactionRule::reactant_container_type::const_iterator j(reactants.begin());
            j != reactants.end(); j++, cnt++)
        {
            r[cnt] = static_cast<double>(world_->get_value_exact(*j));
        }

        cnt = 0;

        ODERatelaw::state_container_type p(products.size());
        for(ODEReactionRule::reactant_container_type::const_iterator j(products.begin());
            j != products.end(); j++, cnt++)
        {
            p[cnt] = static_cast<double>(world_->get_value_exact(*j));
        }

        return rr.get_ratelaw()->deriv_func(r, p, world_->volume(), world_->t(), rr);
    }

protected:
    std::pair<deriv_func, jacobi_func> generate_system() const;
protected:
    // boost::shared_ptr<ODENetworkModel> model_;
    // boost::shared_ptr<ODEWorld_New> world_;
    Real dt_;
    // Integer num_steps_;
    Real abs_tol_, rel_tol_;
    ODESolverType solver_type_;
};

} // ode

} // ecell4

//XXX: #include <cstring>
//XXX: #include <vector>
//XXX: #include <numeric>
//XXX: #include <map>
//XXX: #include <boost/shared_ptr.hpp>
//XXX: #include <boost/scoped_ptr.hpp>
//XXX: #include <boost/scoped_array.hpp>
//XXX: 
//XXX: #include <boost/numeric/ublas/vector.hpp>
//XXX: #include <boost/numeric/ublas/matrix.hpp>
//XXX: 
//XXX: #include <ecell4/core/exceptions.hpp>
//XXX: #include <ecell4/core/types.hpp>
//XXX: #include <ecell4/core/get_mapper_mf.hpp>
//XXX: 
//XXX: #include <ecell4/core/SimulatorBase.hpp>
//XXX: #include <ecell4/core/NetworkModel.hpp>
//XXX: //#include "ODEWorld.hpp"
//XXX: #include <ecell4/core/Species.hpp>
//XXX: #include <ecell4/core/Context.hpp>
//XXX: #include <ecell4/core/Real3.hpp>
//XXX: #include <ecell4/core/Space.hpp>
//XXX: #include <ecell4/core/Model.hpp>
//XXX: #include <ecell4/core/NetfreeModel.hpp>
//XXX: 
//XXX: #include "ODESimulator.hpp"
//XXX: 
//XXX: namespace ecell4
//XXX: {
//XXX: namespace ode
//XXX: {
//XXX: 
//XXX: class ODEWorld_New
//XXX:     : public Space
//XXX: {
//XXX: public:
//XXX:     typedef std::vector<Real> num_molecules_container_type;
//XXX:     typedef std::vector<Species> species_container_type;
//XXX:     typedef utils::get_mapper_mf<
//XXX:         Species, num_molecules_container_type::size_type>::type species_map_type;
//XXX: 
//XXX: public:
//XXX:     ODEWorld_New(const Real3& edge_lengths = Real3(1, 1, 1))
//XXX:         : t_(0.0)
//XXX:     {
//XXX:         reset(edge_lengths);
//XXX:     }
//XXX: 
//XXX:     ODEWorld_New(const std::string& filename)
//XXX:         : t_(0.0)
//XXX:     {
//XXX:         reset(Real3(1, 1, 1));
//XXX:         //this->load(filename);
//XXX:     }
//XXX: 
//XXX:     // Space Traints
//XXX:     const Real t() const
//XXX:     {
//XXX:         return t_;
//XXX:     }
//XXX: 
//XXX:     void set_t(const Real& t)
//XXX:     {
//XXX:         if (t < 0.0)
//XXX:         {
//XXX:             throw std::invalid_argument("the time must be positive.");
//XXX:         }
//XXX:         t_ = t;
//XXX:     }
//XXX: 
//XXX:     const Real3& edge_lengths() const
//XXX:     {
//XXX:         return edge_lengths_;
//XXX:     }
//XXX: 
//XXX:     void reset(const Real3& edge_lengths)
//XXX:     {
//XXX:         t_ = 0.0;
//XXX:         index_map_.clear();
//XXX:         num_molecules_.clear();
//XXX:         species_.clear();
//XXX: 
//XXX:         for (Real3::size_type dim(0); dim < 3; ++dim)
//XXX:         {
//XXX:             if (edge_lengths[dim] <= 0)
//XXX:             {
//XXX:                 throw std::invalid_argument("the edge length must be positive.");
//XXX:             }
//XXX:         }
//XXX: 
//XXX:         edge_lengths_ = edge_lengths;
//XXX:         volume_ = edge_lengths[0] * edge_lengths[1] * edge_lengths[2];
//XXX:     }
//XXX: 
//XXX:     const Real volume() const
//XXX:     {
//XXX:         return volume_;
//XXX:     }
//XXX: 
//XXX:     void set_volume(const Real& volume)
//XXX:     {
//XXX:         if (volume <= 0.0)
//XXX:         {
//XXX:             throw std::invalid_argument("The volume must be positive.");
//XXX:         }
//XXX: 
//XXX:         volume_ = volume;
//XXX:         const Real L(cbrt(volume));
//XXX:         edge_lengths_ = Real3(L, L, L);
//XXX:     }
//XXX: 
//XXX:     Integer num_molecules(const Species& sp) const
//XXX:     {
//XXX:         return static_cast<Integer>(get_value(sp));
//XXX:     }
//XXX: 
//XXX:     Integer num_molecules_exact(const Species& sp) const
//XXX:     {
//XXX:         return static_cast<Integer>(get_value_exact(sp));
//XXX:     }
//XXX: 
//XXX:     std::vector<Species> list_species() const
//XXX:     {
//XXX:         return species_;
//XXX:     }
//XXX: 
//XXX:     // CompartmentSpace member functions
//XXX: 
//XXX:     void add_molecules(const Species& sp, const Real& num)
//XXX:     {
//XXX:         species_map_type::const_iterator i(index_map_.find(sp));
//XXX:         if (i == index_map_.end())
//XXX:         {
//XXX:             reserve_species(sp);
//XXX:             i = index_map_.find(sp);
//XXX:         }
//XXX: 
//XXX:         num_molecules_[(*i).second] += num;
//XXX:     }
//XXX: 
//XXX:     void remove_molecules(const Species& sp, const Real& num)
//XXX:     {
//XXX:         species_map_type::const_iterator i(index_map_.find(sp));
//XXX:         if (i == index_map_.end())
//XXX:         {
//XXX:             throw NotFound("Species not found");
//XXX:         }
//XXX: 
//XXX:         num_molecules_[(*i).second] -= num;
//XXX:     }
//XXX: 
//XXX:     // Optional members
//XXX: 
//XXX:     Real get_value(const Species& sp) const
//XXX:     {
//XXX:         SpeciesExpressionMatcher sexp(sp);
//XXX:         Real retval(0);
//XXX:         for (species_map_type::const_iterator i(index_map_.begin());
//XXX:             i != index_map_.end(); ++i)
//XXX:         {
//XXX:             if (sexp.match((*i).first))
//XXX:             {
//XXX:                 do
//XXX:                 {
//XXX:                     retval += num_molecules_[(*i).second];
//XXX:                 } while (sexp.next());
//XXX:             }
//XXX:         }
//XXX:         return retval;
//XXX:     }
//XXX: 
//XXX:     Real get_value_exact(const Species& sp) const
//XXX:     {
//XXX:         species_map_type::const_iterator i(index_map_.find(sp));
//XXX:         if (i == index_map_.end())
//XXX:         {
//XXX:             // throw NotFound("Species not found");
//XXX:             return 0.0;
//XXX:         }
//XXX: 
//XXX:         return num_molecules_[(*i).second];
//XXX:     }
//XXX: 
//XXX:     void set_value(const Species& sp, const Real& num)
//XXX:     {
//XXX:         species_map_type::const_iterator i(index_map_.find(sp));
//XXX:         if (i == index_map_.end())
//XXX:         {
//XXX:             reserve_species(sp);
//XXX:             i = index_map_.find(sp);
//XXX:         }
//XXX: 
//XXX:         num_molecules_[(*i).second] = num;
//XXX:     }
//XXX: 
//XXX:     //FIXME Following 2 member functions are not implemented for now!
//XXX:     void save(const std::string& filename) const;
//XXX:     void load(const std::string& filename);
//XXX: 
//XXX:     bool has_species(const Species& sp) const
//XXX:     {
//XXX:         species_map_type::const_iterator i(index_map_.find(sp));
//XXX:         return (i != index_map_.end());
//XXX:     }
//XXX: 
//XXX:     void reserve_species(const Species& sp)
//XXX:     {
//XXX:         species_map_type::const_iterator i(index_map_.find(sp));
//XXX:         if (i != index_map_.end())
//XXX:         {
//XXX:             throw AlreadyExists("Species already exists");
//XXX:         }
//XXX: 
//XXX:         index_map_.insert(std::make_pair(sp, num_molecules_.size()));
//XXX:         species_.push_back(sp);
//XXX:         num_molecules_.push_back(0);
//XXX:     }
//XXX: 
//XXX:     void release_species(const Species& sp)
//XXX:     {
//XXX:         species_map_type::iterator i(index_map_.find(sp));
//XXX:         if (i == index_map_.end())
//XXX:         {
//XXX:             throw NotFound("Species not found");
//XXX:         }
//XXX: 
//XXX:         species_map_type::mapped_type
//XXX:             idx((*i).second), last_idx(num_molecules_.size() - 1);
//XXX:         if (idx != last_idx)
//XXX:         {
//XXX:             const species_container_type::size_type
//XXX:                 idx_(static_cast<species_container_type::size_type>(idx)),
//XXX:                 last_idx_(
//XXX:                     static_cast<species_container_type::size_type>(last_idx));
//XXX:             const Species& last_sp(species_[last_idx_]);
//XXX:             species_[idx_] = last_sp;
//XXX:             num_molecules_[idx] = num_molecules_[last_idx];
//XXX:             index_map_[last_sp] = idx;
//XXX:         }
//XXX: 
//XXX:         species_.pop_back();
//XXX:         num_molecules_.pop_back();
//XXX:         index_map_.erase(sp);
//XXX:     }
//XXX: 
//XXX:     void bind_to(boost::shared_ptr<Model> model);
//XXX:     //void bind_to(boost::shared_ptr<NetworkModel> model);
//XXX: 
//XXX:     //boost::shared_ptr<NetworkModel> lock_model() const
//XXX:     boost::shared_ptr<Model> lock_model() const
//XXX:     {
//XXX:         //if (generated_)
//XXX:         //{
//XXX:         //    return generated_;
//XXX:         //}
//XXX:         //else
//XXX:         //{
//XXX:         //    return model_.lock();
//XXX:         //}
//XXX:         return model_.lock();
//XXX:     }
//XXX:     //
//XXX:     void add_molecules(const Species& sp, const Integer& num,
//XXX:         const boost::shared_ptr<Shape> shape)
//XXX:     {
//XXX:         add_molecules(sp, num);
//XXX:     }
//XXX: 
//XXX:     std::pair<std::pair<ParticleID, Particle>, bool> new_particle(const Particle& p)
//XXX:     {
//XXX:         add_molecules(p.species(), 1);
//XXX:         return std::make_pair(std::make_pair(ParticleID(), p), true);
//XXX:     }
//XXX: 
//XXX:     std::pair<std::pair<ParticleID, Particle>, bool> new_particle(
//XXX:         const Species& sp, const Real3& pos)
//XXX:     {
//XXX:         add_molecules(sp, 1);
//XXX:         return std::make_pair(
//XXX:             std::make_pair(ParticleID(), Particle(sp, pos, 0.0, 0.0)), true);
//XXX:     }
//XXX: 
//XXX:     //Real evaluate(const ReactionRule& rr) const
//XXX:     //{
//XXX:     //    ODEReactionRule oderr(rr);
//XXX:     //    return evaluate(oderr);
//XXX:     //}
//XXX: 
//XXX:     //Real evaluate(ODEReactionRule& rr) const;
//XXX: 
//XXX: protected:
//XXX:     Real3 edge_lengths_;
//XXX:     Real volume_;
//XXX:     Real t_;
//XXX: 
//XXX:     num_molecules_container_type num_molecules_;
//XXX:     species_container_type species_;
//XXX:     species_map_type index_map_;
//XXX: 
//XXX:     boost::weak_ptr<Model> model_;
//XXX:     //boost::shared_ptr<NetworkModel> generated_;
//XXX: };
//XXX: 
//XXX: class ODESimulator_New
//XXX:     : public SimulatorBase<Model, ODEWorld_New>
//XXX: {
//XXX: public:
//XXX:     typedef SimulatorBase<Model, ODEWorld_New> base_type;
//XXX: public:
//XXX:     typedef boost::numeric::ublas::vector<Real> state_type;
//XXX:     typedef boost::numeric::ublas::matrix<Real> matrix_type;
//XXX:     typedef std::vector<state_type::size_type> index_container_type;
//XXX:     typedef std::vector<Real> coefficient_container_type;
//XXX:     typedef ReactionRule reaction_container_type;
//XXX:     typedef Model::reaction_rule_container_type reaction_rule_container_type;
//XXX: 
//XXX:     struct mapped_reaction_type
//XXX:     {
//XXX:         index_container_type reactants;
//XXX:         coefficient_container_type reactant_coefficients;
//XXX:         index_container_type products;
//XXX:         coefficient_container_type product_coefficients;
//XXX:         Real k;
//XXX:         //boost::weak_ptr<ODERatelaw> ratelaw;
//XXX:         const ReactionRule *raw;
//XXX:         boost::weak_ptr<ReactionRuleDescriptor> rr_descriptor;
//XXX:     };
//XXX:     typedef std::vector<mapped_reaction_type> mapped_reaction_container_type;
//XXX: 
//XXX:     class deriv_func 
//XXX:     {
//XXX:     public:
//XXX:         deriv_func(const mapped_reaction_container_type &reactions, const Real &volume)
//XXX:             : reactions_(reactions), volume_(volume), vinv_(1./volume)
//XXX:         {;}
//XXX:         void operator()(const state_type &x, state_type &dxdt, const double &t)
//XXX:         {
//XXX:             std::fill(dxdt.begin(), dxdt.end(), 0.);
//XXX:             for(mapped_reaction_container_type::const_iterator mapped_rr_it(reactions_.begin());
//XXX:                     mapped_rr_it != reactions_.end(); mapped_rr_it++) 
//XXX:             {   // Calculate the flux of each ReactionRule
//XXX:                 //  FIXME.  Currently, this can handle only the mass action
//XXX:                 Real flux = (mapped_rr_it->k * volume_);
//XXX:                 //for (index_container_type::const_iterator r = mapped_rr_it->reactants.begin(); 
//XXX:                 //        r != mapped_rr_it->reactants.end(); r++) {
//XXX:                 //    //flux *= std::pow( x[*r] * vinv_, 1. );
//XXX:                 //    flux *=  x[*r] * vinv_;
//XXX:                 //}
//XXX:                 for( index_container_type::size_type i = 0; i < mapped_rr_it->reactants.size(); i++) {
//XXX:                     index_container_type::value_type reactant_index = mapped_rr_it->reactants[i];
//XXX:                     flux *= std::pow(x[reactant_index] * vinv_, mapped_rr_it->reactant_coefficients[i]);
//XXX:                 }
//XXX: 
//XXX:                 // Merge each reactions's flux 
//XXX:                 //for(index_container_type::const_iterator r = mapped_rr_it->reactants.begin();
//XXX:                 //        r != mapped_rr_it->reactants.end(); r++) 
//XXX:                 //{
//XXX:                 //    dxdt[*r] -= (flux); // 
//XXX:                 //    //dxdt[*j] -= (flux * (double)i->reactant_coefficients[nth]);
//XXX:                 //}
//XXX:                 //for(index_container_type::const_iterator p = mapped_rr_it->products.begin();
//XXX:                 //        p != mapped_rr_it->products.end(); p++)
//XXX:                 //{
//XXX:                 //    dxdt[*p] += (flux);
//XXX:                 //}
//XXX:                 for( index_container_type::size_type i = 0; i < mapped_rr_it->reactants.size(); i++) {
//XXX:                     index_container_type::value_type reactant_index = mapped_rr_it->reactants[i];
//XXX:                     dxdt[reactant_index] -= (flux * (double)mapped_rr_it->reactant_coefficients[i]);
//XXX:                 }
//XXX:                 for( index_container_type::size_type i = 0; i < mapped_rr_it->products.size(); i++) {
//XXX:                     index_container_type::value_type product_index = mapped_rr_it->products[i];
//XXX:                     dxdt[product_index] += (flux * (double)mapped_rr_it->product_coefficients[i]);
//XXX:                 }
//XXX:             }
//XXX:             //std::cout << "========================================" << std::endl;
//XXX:             //std::cout << "dxdt dumping at t = " << t << std::endl;
//XXX:             //for (state_type::const_iterator it = dxdt.begin();
//XXX:             //        it != dxdt.end(); it++) 
//XXX:             //{
//XXX:             //    std::cout << *it << std::endl;
//XXX:             //}
//XXX:             return;
//XXX:         }
//XXX:     protected:
//XXX:         const mapped_reaction_container_type reactions_;
//XXX:         const Real volume_;
//XXX:         const Real vinv_;
//XXX:     };
//XXX: 
//XXX:     class jacobi_func
//XXX:     {
//XXX:     public:
//XXX:         jacobi_func(
//XXX:             //const reaction_container_type &reactions, const Real& volume,
//XXX:             const mapped_reaction_container_type &reactions, const Real& volume,
//XXX:             const Real& abs_tol, const Real& rel_tol)
//XXX:             : reactions_(reactions), volume_(volume), vinv_(1.0 / volume), abs_tol_(abs_tol), rel_tol_(rel_tol)
//XXX:         {
//XXX:             ;
//XXX:         }
//XXX:         void operator()(
//XXX:                 const state_type& x, matrix_type& jacobi, const double &t, state_type &dfdt) const
//XXX:         {
//XXX:             //fill 0 into jacobi and dfdt
//XXX:             std::fill(dfdt.begin(), dfdt.end(), 0.0);
//XXX:             std::fill(jacobi.data().begin(), jacobi.data().end(), 0.0);
//XXX: 
//XXX:             // const Real ETA(2.2204460492503131e-16);
//XXX:             const Real SQRTETA(1.4901161193847656e-08);
//XXX:             const Real r0(1.0);
//XXX:             // XXX Following comment-outed code is derived from previous ODESimulator.hpp
//XXX:             //          Should be erased??????
//XXX:             // Real fac(0.0);
//XXX:             // for (std::size_t k(0); k < dfdt.size(); ++k)
//XXX:             // {
//XXX:             //     const Real ewtk(atol + rtol * x[k]);
//XXX:             //     fac = std::max(fac, dfdt[k] * ewtk);
//XXX:             // }
//XXX:             // const Real r0(1000.0 * h * ETA * dfdt.size() * fac);  //XXX: h means the step interval
//XXX:             // {
//XXX:             //     const Real ewtj(atol + rtol * x[j]);
//XXX:             //     const Real dyj(std::max(SQRTETA * abs(x[j]), r0 * ewtj));
//XXX:             // }
//XXX: 
//XXX:             // const Real h(1.0e-8);
//XXX:             // const Real ht(1.0e-10);
//XXX:             const Real ht(1.0e-10);
//XXX: 
//XXX:             //for(mapped_reaction_container_type::const_iterator mapped_rr_it = (reactions_.begin());
//XXX:             //        mapped_rr_it != reactions_.end(); mapped_rr_it++) 
//XXX:             //{
//XXX:             //    // 1. Calculate the current flux. This is the same as the flux calculated in deriv_func.
//XXX:             //    Real flux0 = (mapped_rr_it->k * volume_);
//XXX:             //    for (index_container_type::const_iterator r = mapped_rr_it->reactants.begin(); 
//XXX:             //            r != mapped_rr_it->reactants.end(); r++)
//XXX:             //    {
//XXX:             //        //flux *= std::pow( x[*r] * vinv_, 1. );
//XXX:             //        flux0 *=  x[*r] * vinv_;
//XXX:             //    }
//XXX:             //    // 2. Differentiate by time
//XXX:             //    {
//XXX:             //        // In the mass reactio, this become ZERO.
//XXX:             //        //  Currently, not implemented.
//XXX:             //    }
//XXX:             //    // Differentiate by each reactants
//XXX:             //    {
//XXX:             //        for(std::size_t i(0); i < mapped_rr_it->reactants.size(); i++) {
//XXX: 
//XXX:             //        }
//XXX:             //            const Real ewt = abs_tol_ + rel_tol_ * abs(reactants_states[j]);
//XXX:             //            const Real h = std::max(SQRTETA * abs(reactants_states[j]), r0 * ewt);
//XXX:             //    }
//XXX:             //    // Differentiate by each products
//XXX:             //}
//XXX:         }
//XXX:     protected:
//XXX:         const mapped_reaction_container_type reactions_;
//XXX:         const Real volume_;
//XXX:         const Real vinv_;
//XXX:         const Real abs_tol_, rel_tol_;
//XXX:     };
//XXX: 
//XXX:     struct StateAndTimeBackInserter
//XXX:     {
//XXX:         typedef std::vector<state_type> state_container_type;
//XXX:         typedef std::vector<Real> time_container_type;
//XXX:         // Constructor
//XXX:         StateAndTimeBackInserter(state_container_type &states, time_container_type &times)
//XXX:             : m_states(states), m_times(times)
//XXX:         {;}
//XXX: 
//XXX:         void operator()(const state_type &x, double t)
//XXX:         {
//XXX:             m_states.push_back(x);
//XXX:             m_times.push_back(t);
//XXX:         }
//XXX:         state_container_type &m_states;
//XXX:         time_container_type &m_times;
//XXX:     };
//XXX: 
//XXX:     ODESimulator_New(
//XXX:         const boost::shared_ptr<Model>& model,
//XXX:         const boost::shared_ptr<ODEWorld_New> &world,
//XXX:         const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
//XXX:         :base_type(model, world), dt_(inf), 
//XXX:          abs_tol_(1e-6), rel_tol_(1e-6), solver_type_(solver_type)
//XXX:     {;}
//XXX:     ODESimulator_New(
//XXX:         const boost::shared_ptr<ODEWorld_New>& world,
//XXX:         const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
//XXX:         : base_type(world), dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6),
//XXX:           solver_type_(solver_type)
//XXX:     {
//XXX:         initialize();
//XXX:     }
//XXX: 
//XXX:     //ODESimulator_New(
//XXX:     //    const boost::shared_ptr<Model>& model,
//XXX:     //    const boost::shared_ptr<ODEWorld_New>& world,
//XXX:     //    const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
//XXX:     //    : base_type(boost::shared_ptr<NetworkModel>(new NetworkModel(model)), world),
//XXX:     //      dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6), solver_type_(solver_type)
//XXX:     //{
//XXX:     //    initialize();
//XXX:     //}
//XXX: 
//XXX:     void initialize(void)
//XXX:     {
//XXX:         const std::vector<Species> species(model_->list_species());
//XXX:         for(std::vector<Species>::const_iterator it = species.begin(); it != species.end(); it++)
//XXX:         {
//XXX:             if (world_->has_species(*it) != true) {
//XXX:                 world_->reserve_species(*it);
//XXX:             }
//XXX:         }
//XXX:     }
//XXX: 
//XXX:     void step(void)
//XXX:     {
//XXX:         step(next_time());
//XXX:         // FIXME Temporary, comment out
//XXX:         //if ( this->model_->has_network_model() )
//XXX:         //{
//XXX:         //    this->model_->update_model();
//XXX:         //}
//XXX:     }
//XXX:     bool step(const Real &upto);
//XXX: 
//XXX:     Real t(void) const
//XXX:     {
//XXX:         return world_->t();
//XXX:     }
//XXX:     void set_t(const Real &t)
//XXX:     {
//XXX:         world_->set_t(t);
//XXX:     }
//XXX:     Real dt(void) const
//XXX:     {
//XXX:         return this->dt_;
//XXX:     }
//XXX:     void set_dt(const Real &dt)
//XXX:     {
//XXX:         if (dt <= 0)
//XXX:         {
//XXX:             throw std::invalid_argument("The step size must be positive.");
//XXX:         }
//XXX:         dt_ = dt;
//XXX:     }
//XXX: 
//XXX:     Real absolute_tolerance() const
//XXX:     {
//XXX:         return abs_tol_;
//XXX:     }
//XXX: 
//XXX:     void set_absolute_tolerance(const Real abs_tol)
//XXX:     {
//XXX:         if (abs_tol < 0)
//XXX:         {
//XXX:             throw std::invalid_argument("A tolerance must be positive or zero.");
//XXX:         }
//XXX:         abs_tol_ = abs_tol;
//XXX:     }
//XXX: 
//XXX:     Real relative_tolerance() const
//XXX:     {
//XXX:         return rel_tol_;
//XXX:     }
//XXX: 
//XXX:     void set_relative_tolerance(const Real rel_tol)
//XXX:     {
//XXX:         if (rel_tol < 0)
//XXX:         {
//XXX:             throw std::invalid_argument("A tolerance must be positive or zero.");
//XXX:         }
//XXX:         rel_tol_ = rel_tol;
//XXX:     }
//XXX: 
//XXX:     //Real evaluate(const ReactionRule& rr) const
//XXX:     //{
//XXX:     //    return evaluate(ODEReactionRule(rr));
//XXX:     //}
//XXX: 
//XXX:     //Real evaluate(const ODEReactionRule& rr) const
//XXX:     //{
//XXX:     //    if (!rr.has_ratelaw())
//XXX:     //    {
//XXX:     //        // std::cout << "No ratelaw was bound. Return zero." << std::endl;
//XXX:     //        return 0.0;
//XXX:     //    }
//XXX: 
//XXX:     //    const ODEReactionRule::reactant_container_type reactants = rr.reactants();
//XXX:     //    const ODEReactionRule::product_container_type products = rr.products();
//XXX: 
//XXX:     //    ODERatelaw::state_container_type::size_type cnt(0);
//XXX: 
//XXX:     //    ODERatelaw::state_container_type r(reactants.size());
//XXX:     //    for(ODEReactionRule::reactant_container_type::const_iterator j(reactants.begin());
//XXX:     //        j != reactants.end(); j++, cnt++)
//XXX:     //    {
//XXX:     //        r[cnt] = static_cast<double>(world_->get_value_exact(*j));
//XXX:     //    }
//XXX: 
//XXX:     //    cnt = 0;
//XXX: 
//XXX:     //    ODERatelaw::state_container_type p(products.size());
//XXX:     //    for(ODEReactionRule::reactant_container_type::const_iterator j(products.begin());
//XXX:     //        j != products.end(); j++, cnt++)
//XXX:     //    {
//XXX:     //        p[cnt] = static_cast<double>(world_->get_value_exact(*j));
//XXX:     //    }
//XXX: 
//XXX:     //    return rr.get_ratelaw()->deriv_func(r, p, world_->volume(), world_->t(), rr);
//XXX:     //}
//XXX: 
//XXX: protected:
//XXX:     //std::pair<deriv_func, jacobi_func> generate_system() const;
//XXX:     deriv_func generate_system() const;
//XXX: 
//XXX: protected:
//XXX:     Real dt_;
//XXX:     Real abs_tol_, rel_tol_;
//XXX:     ODESolverType solver_type_;
//XXX: };
//XXX: 
//XXX: } // ode
//XXX: } // ecell4

#endif  // ECELL4_ODE_ODE_SIMULATOR2_HPP
