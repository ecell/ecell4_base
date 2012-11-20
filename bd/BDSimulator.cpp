#include "BDSimulator.hpp"


namespace ecell4
{

namespace bd
{

void BDSimulator::step()
{
    ++(*state_).num_steps;
}

bool BDSimulator::step(Real const& upto)
{
    Real const t_(t()), dt_(dt());
    Real const next_time(t_ + dt_);

    if (upto > next_time)
    {
        step();
        return false;
    }
    else
    {
        set_dt(next_time - t_);
        step();
        set_dt(dt_);
        return true;
    }
}

} // bd

} // ecell4
