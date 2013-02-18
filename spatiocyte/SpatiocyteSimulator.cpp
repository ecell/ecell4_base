#include "SpatiocyteSimulator.hpp"


namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteSimulator::step()
{
    (*world_).step();

    ++num_steps_;
}

bool SpatiocyteSimulator::step(const Real& upto)
{
    const bool retval((*world_).step(upto));
    ++num_steps_;
    return retval;
}

} // spatiocyte

} // ecell4
