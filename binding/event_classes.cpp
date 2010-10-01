#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "Event.hpp"
#include "PythonEvent.hpp"
#include "EventScheduler.hpp"

namespace binding {

typedef EventScheduler<EGFRDSimulatorTraits::time_type> EventSchedulerImpl;

class PythonEvent: public EventSchedulerImpl::Event
{
public:
    typedef EventSchedulerImpl::Event base_type;
    typedef base_type::time_type time_type;

public:
    virtual ~PythonEvent() {}

    PythonEvent(time_type const& time): base_type(time) {}

    PythonEvent(time_type const& time, boost::python::object const& data)
        : base_type(time), data_(data) {}

    boost::python::object const& data() const
    {
        return data_;
    }

private:
    boost::python::object data_;
};

void register_event_class()
{
    register_event_class<EventSchedulerImpl::Event>("Event");
}

void register_python_event_class()
{
    register_python_event_class<PythonEvent>("PythonEvent");
}

void register_event_scheduler_class()
{
    register_event_scheduler_class<EventSchedulerImpl>("EventScheduler");
}

} //namespace binding
