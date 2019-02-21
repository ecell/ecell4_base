#include <pybind11/pybind11.h>
#include <ecell4/spatiocyte/SpatiocyteWorld.hpp>
#include "python_api.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
using namespace ecell4::python_api;

PYBIND11_MODULE(ecell4, m) {
    py::module m_bd         = m.def_submodule("bd",         "A submodule of ecell4");
    py::module m_egfrd      = m.def_submodule("egfrd",      "A submodule of ecell4");
    py::module m_gillespie  = m.def_submodule("gillespie",  "A submodule of ecell4");
    py::module m_meso       = m.def_submodule("meso",       "A submodule of ecell4");
    py::module m_ode        = m.def_submodule("ode",        "A submodule of ecell4");
    py::module m_spatiocyte = m.def_submodule("spatiocyte", "A submodule of ecell4");

    setup_module(m);
    setup_gillespie_module(m_gillespie);
    setup_spatiocyte_module(m_spatiocyte);
}
