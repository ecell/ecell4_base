#ifndef ECELL4_PYTHON_API_SIMULATOR_FACTORY_HPP
#define ECELL4_PYTHON_API_SIMULATOR_FACTORY_HPP

#include <pybind11/pybind11.h>
#include <ecell4/core/SimulatorFactory.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

template<class Factory>
static inline
void define_factory_functions(py::class_<Factory>& factory)
{
    using world_type = typename Factory::world_type;
    using simulator_type = typename Factory::simulator_type;

    factory
        .def("world",
            (world_type* (Factory::*)(const Real3&) const) &Factory::world,
            py::arg("edge_lengths") = ones())
        .def("world",
            (world_type* (Factory::*)(const Real) const) &Factory::world,
            py::arg("volume"))
        .def("world",
            (world_type* (Factory::*)(const std::string&) const) &Factory::world,
            py::arg("filename"))
        .def("world",
            (world_type* (Factory::*)(const boost::shared_ptr<Model>&) const) &Factory::world,
            py::arg("model"))
        .def("simulator",
            (simulator_type* (Factory::*)(const boost::shared_ptr<world_type>&) const) &Factory::simulator,
            py::arg("world"))
        .def("simulator",
            (simulator_type* (Factory::*)(const boost::shared_ptr<world_type>&, const boost::shared_ptr<Model>&) const) &Factory::simulator,
            py::arg("world"), py::arg("model"));
}

}

}

#endif /* ECELL4_PYTHON_API_SIMULATOR_FACTORY_HPP */
