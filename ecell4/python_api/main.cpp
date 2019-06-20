#include <pybind11/pybind11.h>
#include <ecell4/spatiocyte/SpatiocyteWorld.hpp>
#include "python_api.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
using namespace ecell4::python_api;

PYBIND11_MODULE(ecell4_base, m) {
    py::module m_core       = m.def_submodule("core",       "A submodule of ecell4_base");
    py::module m_bd         = m.def_submodule("bd",         "A submodule of ecell4_base");
    py::module m_egfrd      = m.def_submodule("egfrd",      "A submodule of ecell4_base");
    py::module m_gillespie  = m.def_submodule("gillespie",  "A submodule of ecell4_base");
    py::module m_meso       = m.def_submodule("meso",       "A submodule of ecell4_base");
    py::module m_ode        = m.def_submodule("ode",        "A submodule of ecell4_base");
    py::module m_spatiocyte = m.def_submodule("spatiocyte", "A submodule of ecell4_base");

    setup_module(m_core);
    setup_bd_module(m_bd);
    setup_egfrd_module(m_egfrd);
    setup_gillespie_module(m_gillespie);
    setup_meso_module(m_meso);
    setup_ode_module(m_ode);
    setup_spatiocyte_module(m_spatiocyte);
}
