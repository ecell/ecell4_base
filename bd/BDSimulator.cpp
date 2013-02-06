#include "BDSimulator.hpp"


namespace ecell4
{

namespace bd
{

void BDSimulator::step()
{
    {
        BDPropagator propagator(*model_, *world_, rng(), dt());
        while (propagator())
        {
            ; // do nothing here
        }
    }

    set_t(t() + dt());
    ++num_steps_;
}

bool BDSimulator::step(Real const& upto)
{
    Real const t0(t()), dt0(dt());
    Real const next_time(t0 + dt0);

    if (upto <= t0)
    {
        return false;
    }

    if (upto >= next_time)
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

} // bd

} // ecell4
