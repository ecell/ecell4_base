#ifndef ECELL4_PYTHON_API_OBSERVER_HPP
#define ECELL4_PYTHON_API_OBSERVER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <ecell4/core/observers.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

    template<class Base = ecell4::Observer>
    class PyObserver: public Base
    {
    public:
        using Base::Base;

        const Real next_time() const override
        {
            PYBIND11_OVERLOAD(const Real, Base, next_time,);
        }

        void reset() override
        {
            PYBIND11_OVERLOAD(void, Base, reset,);
        }
    };

    class FixedIntervalPythonHooker
        : public FixedIntervalObserver
    {
    public:
        using base_type = FixedIntervalObserver;
        using callback_t = std::function<bool (const boost::shared_ptr<WorldInterface>&, bool)>;

        FixedIntervalPythonHooker(const Real& dt, callback_t callback)
            : base_type(dt), callback_(callback)
        {
        }

        bool fire(const Simulator* sim, const boost::shared_ptr<WorldInterface>& world) override
        {
            return callback_(world, sim->check_reaction()) && base_type::fire(sim, world);
        }

    protected:
        callback_t callback_;

    };

}

}

#endif /* ECELL4_PYTHON_API_OBSERVER_HPP */
