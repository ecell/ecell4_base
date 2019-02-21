#ifndef ECELL4_PYTHON_API_SIMULATOR_HPP
#define ECELL4_PYTHON_API_SIMULATOR_HPP

#include <pybind11/pybind11.h>
#include <ecell4/core/Simulator.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

template<class Base=Simulator>
class PySimulator : public Base
{
public:
    using Base::Base;

    void initialize() override
    {
        PYBIND11_OVERLOAD_PURE(void, Base, initialize,);
    }

    Real t() const override
    {
        PYBIND11_OVERLOAD_PURE(Real, Base, t,);
    }

    Real dt() const override
    {
        PYBIND11_OVERLOAD_PURE(Real, Base, dt,);
    }

    void set_dt(const Real& dt) override
    {
        PYBIND11_OVERLOAD_PURE(void, Base, set_dt, dt);
    }

    Integer num_steps() const override
    {
        PYBIND11_OVERLOAD_PURE(Integer, Base, num_steps,);
    }

    void step() override
    {
        PYBIND11_OVERLOAD_PURE(void, Base, step,);
    }

    bool step(const Real& upto) override
    {
        PYBIND11_OVERLOAD_PURE(bool, Base, step, upto);
    }

    bool check_reaction() const override
    {
        PYBIND11_OVERLOAD(bool, Base, check_reaction,);
    }
};

template<class Base>
class PySimulatorImpl : public PySimulator<Base>
{
public:

    using PySimulator<Base>::PySimulator;

    void initialize() override
    {
        PYBIND11_OVERLOAD(void, Base, initialize,);
    }

    Real t() const override
    {
        PYBIND11_OVERLOAD(Real, Base, t,);
    }

    Real dt() const override
    {
        PYBIND11_OVERLOAD(Real, Base, dt,);
    }

    void set_dt(const Real& dt) override
    {
        PYBIND11_OVERLOAD(void, Base, set_dt, dt);
    }

    Integer num_steps() const override
    {
        PYBIND11_OVERLOAD(Integer, Base, num_steps,);
    }

    void step() override
    {
        PYBIND11_OVERLOAD(void, Base, step,);
    }

    bool step(const Real& upto) override
    {
        PYBIND11_OVERLOAD(bool, Base, step, upto);
    }
};

template<class S, class... Others>
static inline
void define_simulator_functions(py::class_<S, Others...>& simulator)
{
    simulator
        .def("last_reactions", &S::last_reactions)
        .def("model", &S::model)
        .def("world", &S::world)
        .def("run",
            (void (S::*)(const Real&, const bool)) &S::run,
            py::arg("duration"), py::arg("is_dirty") = true)
        .def("run",
            (void (S::*)(const Real&, const boost::shared_ptr<Observer>&, const bool)) &S::run,
            py::arg("duration"), py::arg("observer"), py::arg("is_dirty") = true)
        .def("run",
            (void (S::*)(const Real&, std::vector<boost::shared_ptr<Observer>>, const bool)) &S::run,
            py::arg("duration"), py::arg("observers"), py::arg("is_dirty") = true)
        ;
}

}

}

#endif /* ECELL4_PYTHON_API_SIMULATOR_HPP */
