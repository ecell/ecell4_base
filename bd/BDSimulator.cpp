#include "BDSimulator.hpp"


namespace ecell4
{

namespace bd
{

void BDSimulator::step()
{
    {
        BDPropagator propagator(*model_, *world_, rng_, dt_);
        while (propagator())
        {
            ; // do nothing here
        }
    }

    set_t(t() + dt_);
    ++num_steps_;
}

bool BDSimulator::step(Real const& upto)
{
    Real const t0(t()), dt0(dt_);
    Real const next_time(t0 + dt0);

    if (upto > next_time)
    {
        step();
        return false;
    }
    else
    {
        set_dt(next_time - t0);
        step();
        set_dt(dt0);
        return true;
    }
}

} // bd

} // ecell4
