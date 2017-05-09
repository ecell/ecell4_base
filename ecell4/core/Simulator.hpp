#ifndef ECELL4_SIMULATOR_HPP
#define ECELL4_SIMULATOR_HPP

#include <vector>
#include "types.hpp"
#include "ReactionRule.hpp"


namespace ecell4
{

class Simulator
{
public:

    virtual ~Simulator()
    {
        ; // do nothing
    }

    // SimulatorTraits

    /**
     * initialize
     */
     virtual void initialize() = 0;

    /**
     * get current time.
     * @return time Real
     */
    virtual Real t() const = 0;

    /**
     * get step interval.
     * @return dt Real
     */
    virtual Real dt() const = 0;

    /**
     * set step interval.
     */
    virtual void set_dt(const Real& dt) = 0;

    /**
     * get the number of steps.
     * @return the number of steps Integer
     */
    virtual Integer num_steps() const = 0;

    /**
     * step.
     */
    virtual void step() = 0;

    /**
     * step and return true if the next time is less than upto.
     * if not, step till upto and return false.
     * @return if the simulator does not rearch upto
     */
    virtual bool step(const Real& upto) = 0;

    /**
     * get next time (t + dt).
     * @return next time Real
     */
    inline Real next_time() const
    {
        return t() + dt();
    }

    /**
     * @return if any reaction occurs at the last step or not
     */
    virtual bool check_reaction() const
    {
        return false;
    }
};

}

#endif /* ECELL4_SIMULATOR_HPP */
