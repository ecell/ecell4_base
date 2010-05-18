#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "ParticleContainer.hpp"
#include "binding_common.hpp"

namespace binding {

void register_particle_container_class()
{
    register_particle_container_class<ParticleContainer>("ParticleContainer");
}

} // namespace binding
