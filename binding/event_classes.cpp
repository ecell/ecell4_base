#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "Event.hpp"
#include "EventScheduler.hpp"

namespace binding {

typedef EventScheduler<EGFRDSimulatorTraits::time_type> EventSchedulerImpl;

void register_event_class()
{
    register_event_class<EventSchedulerImpl::Event>("Event");
}

void register_event_scheduler_class()
{
    register_event_scheduler_class<EventSchedulerImpl>("EventScheduler");
}

} //namespace binding
