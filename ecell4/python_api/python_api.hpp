#ifndef ECELL4_PYTHON_API_PYTHON_API_HPP
#define ECELL4_PYTHON_API_PYTHON_API_HPP

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include "type_caster.hpp"

namespace ecell4
{

namespace python_api
{

void setup_module(pybind11::module& m);
void setup_bd_module(pybind11::module& m);
void setup_spatiocyte_module(pybind11::module& m);
void setup_gillespie_module(pybind11::module& m);

}

}
#endif /* ECELL4_PYTHON_API_PYTHON_API_HPP */
