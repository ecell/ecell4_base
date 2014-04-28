#ifndef __ECELL4_SIMULATOR_HPP
#define __ECELL4_SIMULATOR_HPP

#include "types.hpp"


namespace ecell4
{

class SimulatorBase
{
public:

    virtual ~SimulatorBase()
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
};

template <typename Tmodel_, typename Tworld_>
class Simulator
    : public SimulatorBase
{
public:

    typedef Tmodel_ model_type;
    typedef Tworld_ world_type;

public:

    Simulator(boost::shared_ptr<model_type> model,
        boost::shared_ptr<world_type> world)
        : model_(model), world_(world), num_steps_(0)
    {
        world_->bind_to(model_);
    }

    Simulator(boost::shared_ptr<world_type> world)
        : world_(world), num_steps_(0)
    {
        if (boost::shared_ptr<model_type> bound_model = world_->lock_model())
        {
            model_ = bound_model;
        }
        else
        {
            throw std::invalid_argument("A world must be bound to a model.");
        }
    }

    virtual ~Simulator()
    {
        ; // do nothing
    }

    boost::shared_ptr<model_type> model() const
    {
        return model_;
    }

    boost::shared_ptr<world_type> world() const
    {
        return world_;
    }

    /**
     * get the number of steps.
     * @return the number of steps Integer
     */
    Integer num_steps() const
    {
        return num_steps_;
    }

protected:

    boost::shared_ptr<model_type> model_;
    boost::shared_ptr<world_type> world_;
    Integer num_steps_;
};

}

#endif /* __ECELL4_SIMULATOR_HPP */
