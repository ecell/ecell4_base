#include <boost/scoped_array.hpp>

#include "EGFRDSimulatorWrapper.hpp"


namespace ecell4
{

namespace egfrd
{

void EGFRDSimulatorWrapper::step()
{
    if ((*world_).num_particles() == 0)
    {
        ; // increment time
        return;
    }

    (*sim_).step();
    (*world_).set_t((*sim_).t());
}

bool EGFRDSimulatorWrapper::step(const Real& upto)
{
    if ((*world_).num_particles() == 0)
    {
        ; // increment time
        return true; // this should be false
    }

    const bool retval((*sim_).step(upto));
    (*world_).set_t((*sim_).t());
    return retval;
}

} // egfrd

} // ecell4
