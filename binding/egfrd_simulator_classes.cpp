#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "EGFRDSimulator.hpp"

namespace binding {

void register_egfrd_simulator_classes()
{
    register_egfrd_simulator_class<EGFRDSimulator>("_EGFRDSimulator");
}

} // namespace binding
