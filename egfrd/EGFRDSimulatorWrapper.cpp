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
    bool retval((*sim_).step(upto));
    (*world_).set_t((*sim_).t());
    return retval;
}

} // egfrd

} // ecell4
