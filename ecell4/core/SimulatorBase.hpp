#ifndef ECELL4_SIMULATOR_BASE_HPP
#define ECELL4_SIMULATOR_BASE_HPP

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

    struct ObserverEvent: Event
    {
        ObserverEvent(
            SimulatorBase<model_type, world_type>* sim, Observer* obs, const Real& t)
            : Event(t), sim_(sim), obs_(obs), running_(true)
        {
            time_ = obs_->next_time();
        }

        virtual ~ObserverEvent()
        {
            ;
        }

        virtual void fire()
        {
            const boost::shared_ptr<Space> space = sim_->world();
            running_ = obs_->fire(sim_, space);
            // running_ = obs_->fire(sim_, sim_->world());
            // running_ = obs_->fire(sim_, static_cast<const Space*>(sim_->world().get()));
            time_ = obs_->next_time();
        }

        bool running() const
        {
            return running_;
        }

    protected:

        SimulatorBase<model_type, world_type>* sim_;
        Observer* obs_;
        bool running_;
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

    void run(const Real& duration)
    {
        const Real upto(t() + duration);
        while (step(upto))
        {
            ; // do nothing
        }
    }

    void run(const Real& duration, const boost::shared_ptr<Observer>& observer)
    {
        std::vector<boost::shared_ptr<Observer> > observers;
        observers.push_back(observer);
        run(duration, observers);
    }

    bool fire_observers(
        const std::vector<boost::shared_ptr<Observer> >::iterator begin,
        const std::vector<boost::shared_ptr<Observer> >::iterator end)
    {
        bool retval = true;
        for (std::vector<boost::shared_ptr<Observer> >::iterator
            i(begin); i != end; ++i)
        {
            // if (!(*i)->fire(this, static_cast<const Space*>(world_.get())))
            if (!(*i)->fire(this, world_))
            {
                retval = false;
            }
        }
        return retval;
    }

    void run(const Real& duration, std::vector<boost::shared_ptr<Observer> > observers)
    {
        const Real upto(t() + duration);

        std::vector<boost::shared_ptr<Observer> >::iterator
            offset(std::partition(
                observers.begin(), observers.end(), observer_every()));

        for (std::vector<boost::shared_ptr<Observer> >::iterator i(observers.begin());
            i != observers.end(); ++i)
        {
            // (*i)->initialize(world_.get());
            (*i)->initialize(world_);
        }

        EventScheduler scheduler;
        // for (std::vector<boost::shared_ptr<Observer> >::const_iterator
        //     i(offset); i != observers.end(); ++i)
        for (std::vector<boost::shared_ptr<Observer> >::const_iterator
            i(observers.begin()); i != observers.end(); ++i)
        {
            scheduler.add(boost::shared_ptr<Event>(
                new ObserverEvent(this, (*i).get(), t())));
        }

        while (true)
        {
            bool running = true;
            while (next_time() < std::min(scheduler.next_time(), upto))
            {
                step();

                if (!fire_observers(observers.begin(), offset))
                {
                    running = false;
                    break;
                }
            }

            if (!running)
            {
                break;
            }
            else if (upto >= scheduler.next_time())
            {
                step(scheduler.next_time());
                if (!fire_observers(observers.begin(), offset))
                {
                    running = false;
                }
                EventScheduler::value_type top(scheduler.pop());
                top.second->fire();
                running = (
                    running && static_cast<ObserverEvent*>(top.second.get())->running());
                scheduler.add(top.second);
                if (!running)
                {
                    break;
                }
            }
            else
            {
                step(upto);
                fire_observers(observers.begin(), offset);
                break;
            }
        }

        for (std::vector<boost::shared_ptr<Observer> >::iterator i(observers.begin());
            i != observers.end(); ++i)
        {
            // (*i)->finalize(world_.get());
            (*i)->finalize(world_);
        }
    }

protected:

    boost::shared_ptr<model_type> model_;
    boost::shared_ptr<world_type> world_;
    Integer num_steps_;
};

}

#endif /* ECELL4_SIMULATOR_BASE_HPP */
