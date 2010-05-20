#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "Plane.hpp"

namespace binding {

void register_plane_class()
{
    register_plane_class<Plane>("Plane");
}

} // namespace binding
