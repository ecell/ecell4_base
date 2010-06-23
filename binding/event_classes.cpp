#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "Event.hpp"
#include "EventScheduler.hpp"

namespace binding {

template<typename Ttime_>
class Event
{
public:
    typedef Ttime_ time_type;

public:
    Event(time_type time, boost::python::object const& callback)
        : time_(time), callback_(callback) {}

    time_type const& time() const
    {
        return time_;
    }

    boost::python::object const& callback() const
    {
        return callback_;
    }

    void fire()
    {
        callback_();
    }

    Event() {}

private:
    time_type time_;
    boost::python::object callback_;
};

typedef Event<EGFRDSimulatorTraits::time_type> EventImpl;
typedef EventScheduler<EventImpl> EventSchedulerImpl;

void register_event_class()
{
    register_event_class<EventImpl>("Event");
}

void register_event_scheduler_class()
{
    register_event_scheduler_class<EventSchedulerImpl>("EventScheduler");
}

} //namespace binding
