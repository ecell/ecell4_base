#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "PyEvent.hpp"
#include "PyEventScheduler.hpp"

namespace binding {

void register_py_event_class()
{
    register_py_event_class<PyEvent>("Event");
}

void register_py_event_scheduler_class()
{
    register_py_event_scheduler_class<PyEventScheduler>("EventScheduler");
}

} //namespace binding
