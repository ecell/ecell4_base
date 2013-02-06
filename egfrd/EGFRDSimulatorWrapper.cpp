#include "EGFRDSimulatorWrapper.hpp"


namespace ecell4
{

namespace egfrd
{

void EGFRDSimulatorWrapper::step()
{
    (*sim_).step();
}

bool EGFRDSimulatorWrapper::step(Real const& upto)
{
    (*sim_).step(upto);
}

} // egfrd

} // ecell4
