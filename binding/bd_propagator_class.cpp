#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "BDPropagator.hpp"
#include "binding_common.hpp"

namespace binding {

void register_bd_propagator_class()
{
    register_bd_propagator_class<BDPropagator>("BDPropagator");
}

} // namespace binding
