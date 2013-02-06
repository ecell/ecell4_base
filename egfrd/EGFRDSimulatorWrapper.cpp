#include "EGFRDSimulatorWrapper.hpp"


namespace ecell4
{

namespace egfrd
{

void EGFRDSimulatorWrapper::step()
{
    ;

    set_t(t() + dt());
    ++num_steps_;
}

bool EGFRDSimulatorWrapper::step(Real const& upto)
{
    return false;
}

} // egfrd

} // ecell4
