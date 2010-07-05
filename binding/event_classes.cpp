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
    Event(time_type time, boost::python::object data)
        : time_(time), data_(data) {}

    time_type const& time() const
    {
        return time_;
    }

    boost::python::object const& data() const
    {
        return data_;
    }

    bool operator==(Event const& rhs) const
    {
        return time_ == rhs.time_ && data_ == rhs.data_;
    }

    bool operator!=(Event const& rhs) const
    {
        return !operator==(rhs);
    }

    Event() {}

private:
    time_type time_;
    boost::python::object data_;
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
