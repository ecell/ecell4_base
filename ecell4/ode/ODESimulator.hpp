#ifndef __ECELL4_ODE_ODE_SIMULATOR_HPP
#define __ECELL4_ODE_ODE_SIMULATOR_HPP

#include <cstring>
#include <vector>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/ModelWrapper.hpp>

#include "ODEWorld.hpp"

#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ecell4
{

namespace ode
{

class ODESimulator
    : public SimulatorBase<NetworkModel, ODEWorld>
{
public:

    typedef SimulatorBase<NetworkModel, ODEWorld> base_type;

public:

    typedef boost::numeric::ublas::vector<double> state_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef std::vector<state_type::size_type> index_container_type;

    typedef struct
    {
        index_container_type reactants;
        index_container_type products;
        Real k;
    } reaction_type;

    typedef std::vector<reaction_type> reaction_container_type;

    class deriv_func
    {
    public:

        deriv_func(const reaction_container_type& reactions, const Real& volume)
            : reactions_(reactions), volume_(volume), vinv_(1.0 / volume)
        {
            ;
        }

        void operator()(const state_type& x, state_type& dxdt, const double& t)
        {
            for (state_type::iterator i(dxdt.begin()); i != dxdt.end(); ++i)
            {
                *i = 0.0;
            }

            for (reaction_container_type::const_iterator
                i(reactions_.begin()); i != reactions_.end(); ++i)
            {
                double flux((*i).k * volume_);
                for (index_container_type::const_iterator
                    j((*i).reactants.begin()); j != (*i).reactants.end(); ++j)
                {
                    flux *= x[*j] * vinv_;
                }

                for (index_container_type::const_iterator
                    j((*i).reactants.begin()); j != (*i).reactants.end(); ++j)
                {
                    dxdt[*j] -= flux;
                }

                for (index_container_type::const_iterator
                    j((*i).products.begin()); j != (*i).products.end(); ++j)
                {
                    dxdt[*j] += flux;
                }
            }
        }

    protected:

        const reaction_container_type reactions_;
        const Real volume_;
        const Real vinv_;
    };

    class jacobi_func
    {
    public:

        jacobi_func(const reaction_container_type& reactions, const Real& volume)
            : reactions_(reactions), volume_(volume), vinv_(1.0 / volume)
        {
            ;
        }

        void operator()(
            const state_type& x, matrix_type& jacobi, const double& t, state_type& dfdt) const
        {
            for (state_type::iterator i(dfdt.begin()); i != dfdt.end(); ++i)
            {
                *i = 0.0;
            }

            for (matrix_type::array_type::iterator i(jacobi.data().begin());
                i != jacobi.data().end(); ++i)
            {
                *i = 0.0;
            }

            for (reaction_container_type::const_iterator
                i(reactions_.begin()); i != reactions_.end(); ++i)
            {
                double flux((*i).k * volume_);
                for (index_container_type::const_iterator
                    j((*i).reactants.begin()); j != (*i).reactants.end(); ++j)
                {
                    flux *= x[*j] * vinv_;
                }

                if (flux == 0)
                {
                    continue;
                }

                for (index_container_type::const_iterator
                    j((*i).reactants.begin()); j != (*i).reactants.end(); ++j)
                {
                    const double partial(flux / x[*j]); //XXX: consider squares

                    for (index_container_type::const_iterator
                        k((*i).reactants.begin()); k != (*i).reactants.end(); ++k)
                    {
                        jacobi(*k, *j) -= partial;
                    }

                    for (index_container_type::const_iterator
                        k((*i).products.begin()); k != (*i).products.end(); ++k)
                    {
                        jacobi(*k, *j) += partial;
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

        state_container_type& m_states;
        time_container_type& m_times;

        StateAndTimeBackInserter(
            state_container_type& states, time_container_type& times)
            : m_states(states), m_times(times)
        {
            ;
        }

        void operator()(const state_type&x, double t)
        {
            m_states.push_back(x);
            m_times.push_back(t);
        }
    };

public:

    ODESimulator(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<ODEWorld> world)
        : base_type(model, world), dt_(inf)
    {
        initialize();
    }

    ODESimulator(boost::shared_ptr<ODEWorld> world)
        : base_type(world), dt_(inf)
    {
        initialize();
    }

    void initialize()
    {
        const std::vector<Species> species(model_->list_species());
        for (std::vector<Species>::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            if (!(*world_).has_species(*i))
            {
                (*world_).reserve_species(*i);
            }
        }
    }

    // SimulatorTraits

    Real t(void) const
    {
        return (*world_).t();
    }

    Real dt() const
    {
        return dt_;
    }

    void step(void)
    {
        step(next_time());
    }

    bool step(const Real& upto);

    // Optional members

    void set_t(const Real& t)
    {
        (*world_).set_t(t);
    }

    void set_dt(const Real& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

protected:

    std::pair<deriv_func, jacobi_func> generate_system() const;

protected:

    Real dt_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_SIMULATOR_HPP */
