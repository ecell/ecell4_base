#ifndef __BD_SIMULATOR_STATE_HPP
#define __BD_SIMULATOR_STATE_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>


namespace ecell4
{

namespace bd
{

struct BDSimulatorState
{
    BDSimulatorState(RandomNumberGenerator& r)
        : rng(r), num_steps(0), dt(0)
    {
        ;
    }

    Real dt;
    Integer num_steps;
    RandomNumberGenerator& rng;
};

} // bd

} // ecell4

#endif /* __BD_SIMULATOR_STATE_HPP */
