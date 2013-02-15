#include "SpatiocyteSimulator.hpp"


namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteSimulator::step()
{
    ; // implement here

    set_t(t() + dt());
    ++num_steps_;
}

bool SpatiocyteSimulator::step(const Real& upto)
{
    const Real t0(t()), dt0(dt()), tnext(next_time());

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
        set_dt(upto - t0);
        step();
        set_dt(dt0);
        return false;
    }
}

} // spatiocyte

} // ecell4
