#ifndef __ECELL4_SIMULATOR_BASE_HPP
#define __ECELL4_SIMULATOR_BASE_HPP

#include <time.h>

#include "Simulator.hpp"
#include "EventScheduler.hpp"
#include "observers.hpp"


namespace ecell4
{

template <typename Tmodel_, typename Tworld_>
class SimulatorBase
    : public Simulator
{
public:

    typedef Tmodel_ model_type;
    typedef Tworld_ world_type;

protected:

    struct ObserverEvent:
        EventScheduler::Event
    {
        ObserverEvent(
            SimulatorBase<model_type, world_type>* sim, Observer* obs, const Real& t)
            : EventScheduler::Event(t), sim_(sim), obs_(obs)
        {
            time_ = obs_->next_time();
        }

        virtual ~ObserverEvent()
        {
            ;
        }

        virtual void fire()
        {
            obs_->fire(sim_, sim_->world().get());
            time_ = obs_->next_time();
        }

    protected:

        SimulatorBase<model_type, world_type>* sim_;
        Observer* obs_;
    };

    struct observer_every
    {
        bool operator()(boost::shared_ptr<Observer> const& val) const
        {
            return val->every();
        }
    };

public:

    SimulatorBase(const boost::shared_ptr<model_type>& model,
        const boost::shared_ptr<world_type>& world)
        : model_(model), world_(world), num_steps_(0)
    {
        world_->bind_to(model_);
    }

    SimulatorBase(const boost::shared_ptr<world_type>& world)
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

    virtual ~SimulatorBase()
    {
        ; // do nothing
    }

    const boost::shared_ptr<model_type>& model() const
    {
        return model_;
    }

    const boost::shared_ptr<world_type>& world() const
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

    virtual Real t() const
    {
        return (*world_).t();
    }

    virtual void set_t(const Real& t)
    {
        (*world_).set_t(t);
    }

    /**
     * set step interval.
     */
    virtual void set_dt(const Real& dt)
    {
        std::cerr << "WARN: set_dt(const Real&) was just ignored." << std::endl;
    }

    Real run(const Real& duration)
    {
        time_t t_start, t_end;
        time(&t_start);

        const Real upto(t() + duration);
        while (step(upto))
        {
            ; // do nothing
        }

        time(&t_end);
        return difftime(t_end, t_start);
    }

    Real run(const Real& duration, const boost::shared_ptr<Observer>& observer)
    {
        std::vector<boost::shared_ptr<Observer> > observers;
        observers.push_back(observer);
        return run(duration, observers);
    }

    Real run(const Real& duration, std::vector<boost::shared_ptr<Observer> > observers)
    {
        time_t t_start, t_end;
        time(&t_start);

        const Real upto(t() + duration);

        std::vector<boost::shared_ptr<Observer> >::iterator
            offset(std::partition(
                observers.begin(), observers.end(), observer_every()));

        for (std::vector<boost::shared_ptr<Observer> >::iterator i(observers.begin());
            i != observers.end(); ++i)
        {
            (*i)->initialize(world_.get());
        }

        EventScheduler scheduler;
        for (std::vector<boost::shared_ptr<Observer> >::const_iterator
            i(offset); i != observers.end(); ++i)
        {
            scheduler.add(boost::shared_ptr<EventScheduler::Event>(
                new ObserverEvent(this, (*i).get(), t())));
        }

        while (true)
        {
            while (next_time() < std::min(upto, scheduler.next_time()))
            {
                step();
                for (std::vector<boost::shared_ptr<Observer> >::iterator
                    i(observers.begin()); i != offset; ++i)
                {
                    (*i)->fire(this, world_.get());
                }
            }

            if (upto >= scheduler.next_time())
            {
                step(scheduler.next_time());
                for (std::vector<boost::shared_ptr<Observer> >::iterator
                    i(observers.begin()); i != offset; ++i)
                {
                    (*i)->fire(this, world_.get());
                }

                EventScheduler::value_type top(scheduler.pop());
                top.second->fire();
                scheduler.add(top.second);
            }
            else
            {
                step(upto);
                for (std::vector<boost::shared_ptr<Observer> >::iterator
                    i(observers.begin()); i != offset; ++i)
                {
                    (*i)->fire(this, world_.get());
                }

                break;
            }
        }

        for (std::vector<boost::shared_ptr<Observer> >::iterator i(observers.begin());
            i != observers.end(); ++i)
        {
            (*i)->finalize(world_.get());
        }

        time(&t_end);
        return difftime(t_end, t_start);
    }

protected:

    boost::shared_ptr<model_type> model_;
    boost::shared_ptr<world_type> world_;
    Integer num_steps_;
};

}

#endif /* __ECELL4_SIMULATOR_BASE_HPP */
