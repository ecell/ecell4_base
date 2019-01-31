#include <pybind11/pybind11.h>
#include <ecell4/spatiocyte/SpatiocyteWorld.hpp>
#include "type_caster.hpp"
#include "core.hpp"

namespace py = pybind11;

namespace {
    using namespace ecell4::spatiocyte;

    void setup_spatiocyte_module(py::module& m)
    {
        py::class_<SpatiocyteWorld>(m, "SpatiocyteWorld")
            .def(py::init<>())
            .def("t", &SpatiocyteWorld::t)
            .def("set_t", &SpatiocyteWorld::set_t)
            .def("save", &SpatiocyteWorld::save)
            .def("load", &SpatiocyteWorld::load)
            .def("volume", &SpatiocyteWorld::volume)
            .def("has_species", &SpatiocyteWorld::has_species);
    }

}

PYBIND11_MODULE(ecell4, m) {
    py::module m_bd         = m.def_submodule("bd",         "A submodule of ecell4");
    py::module m_egfrd      = m.def_submodule("egfrd",      "A submodule of ecell4");
    py::module m_gillespie  = m.def_submodule("gillespie",  "A submodule of ecell4");
    py::module m_meso       = m.def_submodule("meso",       "A submodule of ecell4");
    py::module m_ode        = m.def_submodule("ode",        "A submodule of ecell4");
    py::module m_spatiocyte = m.def_submodule("spatiocyte", "A submodule of ecell4");

    setup_module(m);
    setup_spatiocyte_module(m_spatiocyte);
}
