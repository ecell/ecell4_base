#include "EGFRDSimulatorWrapper.hpp"


namespace ecell4
{

namespace egfrd
{

void EGFRDSimulatorWrapper::step()
{
    (*sim_).step();
    (*world_).set_t((*sim_).t());
}

bool EGFRDSimulatorWrapper::step(Real const& upto)
{
    (*sim_).step(upto);
    (*world_).set_t((*sim_).t());
}

} // egfrd

} // ecell4
