#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Structure.hpp"
#include "ParticleSimulationStructure.hpp"
#include "Surface.hpp"
#include "Region.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "SphericalSurface.hpp"
#include "CuboidalRegion.hpp"
#include "binding_common.hpp"

namespace binding {

void register_structure_classes()
{
    register_structure_class<Structure>("Structure");
    register_particle_simulation_structure_class<ParticleSimulationStructure>("ParticleSimulationStructure");
    register_surface_class<Surface>("Surface");
    register_region_class<Region>("Region");
    register_cuboidal_region_class<CuboidalRegion>("CuboidalRegion");
    register_planar_surface_class<PlanarSurface>("PlanarSurface");
    register_spherical_surface_class<SphericalSurface>("SphericalSurface");
    register_cylindrical_surface_class<CylindricalSurface>("CylindricalSurface");
}

} // namespace binding

