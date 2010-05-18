#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Particle.hpp"
#include "binding_common.hpp"

namespace binding {

void register_particle_class()
{
    ParticleWrapper<Particle>::__register_class("Particle");
}

} // namespace binding
