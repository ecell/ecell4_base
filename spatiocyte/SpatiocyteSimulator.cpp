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
    const Real t0(t()), tnext(next_time());
    // const Real dt0(dt());

    if (upto <= t0)
    {
        return false;
    }

    if (upto >= tnext)
    {
        step();
        return true;
    }
    else
    {
        // set_dt(upto - t0);
        step();
        // set_dt(dt0);
        return false;
    }
}

} // spatiocyte

} // ecell4
