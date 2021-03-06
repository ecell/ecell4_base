#ifndef ECELL4_ODE_ODE_FACTORY2_HPP
#define ECELL4_ODE_ODE_FACTORY2_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include "ODEWorld.hpp"
#include "ODESimulator.hpp"


namespace ecell4
{

namespace ode
{

class ODEFactory:
    public SimulatorFactory<ODEWorld, ODESimulator>
{
public:

    typedef SimulatorFactory<ODEWorld, ODESimulator> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::simulator_type simulator_type;
    typedef ODEFactory this_type;

public:

    ODEFactory(const ODESolverType solver_type = default_solver_type(), const Real dt = default_dt(),
               const Real abs_tol = default_abs_tol(), const Real rel_tol = default_rel_tol())
        : base_type(), solver_type_(solver_type), dt_(dt), abs_tol_(abs_tol), rel_tol_(rel_tol)
    {
        ; // do nothing
    }

    virtual ~ODEFactory()
    {
        ; // do nothing
    }

    static inline const ODESolverType default_solver_type()
    {
        return ROSENBROCK4_CONTROLLER;
    }

    static inline const Real default_dt()
    {
        return std::numeric_limits<Real>::infinity();
    }

    static inline const Real default_abs_tol()
    {
        return 0.0;
    }

    static inline const Real default_rel_tol()
    {
        return 0.0;
    }

    ODEFactory& rng(const std::shared_ptr<RandomNumberGenerator>& rng)
    {
        return (*this);  //XXX: Just for the compatibility
    }

    inline ODEFactory* rng_ptr(const std::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
    }

protected:

    virtual simulator_type* create_simulator(
        const std::shared_ptr<world_type>& w, const std::shared_ptr<Model>& m) const
    {
        simulator_type* sim = new simulator_type(w, m, solver_type_);
        sim->set_dt(dt_);

        if (abs_tol_ > 0)
        {
            sim->set_absolute_tolerance(abs_tol_);
        }

        if (rel_tol_ > 0)
        {
            sim->set_relative_tolerance(rel_tol_);
        }
        return sim;
    }

protected:

    ODESolverType solver_type_;
    Real dt_, abs_tol_, rel_tol_;
};

} // ode

} // ecell4

#endif /* ECELL4_ODE_ODE_FACTORY2_HPP */
