#ifndef __SIMULATOR_HPP
#define __SIMULATOR_HPP

#include "types.hpp"


namespace ecell4
{

class Simulator
{
public:

    virtual Real t() const = 0;
    virtual Integer num_steps() const = 0;

    virtual void step() = 0;
    virtual bool step(Real const& upto) = 0;
};

}

#endif /* __SIMULATOR_HPP */
