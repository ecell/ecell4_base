#ifndef ECELL4_PYTHON_API_SIMULATOR_FACTORY_HPP
#define ECELL4_PYTHON_API_SIMULATOR_FACTORY_HPP

#include <pybind11/pybind11.h>

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
            [](const Factory& self, const Real3& edge_lengths)
            {
                return boost::shared_ptr<world_type>(self.world(edge_lengths));
            },
            py::arg("edge_lengths") = ones())
        .def("world",
            [](const Factory& self, const Real volume)
            {
                return boost::shared_ptr<world_type>(self.world(volume));
            },
            py::arg("volume"))
        .def("world",
            [](const Factory& self, const std::string& filename)
            {
                return boost::shared_ptr<world_type>(self.world(filename));
            },
            py::arg("filename"))
        .def("world",
            [](const Factory& self, const boost::shared_ptr<Model>& model)
            {
                return boost::shared_ptr<world_type>(self.world(model));
            },
            py::arg("model"))
        .def("simulator",
            [](const Factory& self, const boost::shared_ptr<world_type>& world)
            {
                return boost::shared_ptr<simulator_type>(self.simulator(world));
            },
            py::arg("world"))
        .def("simulator",
            [](const Factory& self, const boost::shared_ptr<world_type>& world, const boost::shared_ptr<Model>& model)
            {
                return boost::shared_ptr<simulator_type>(self.simulator(world, model));
            },
            py::arg("world"), py::arg("model"));
}

}

}

#endif /* ECELL4_PYTHON_API_SIMULATOR_FACTORY_HPP */
