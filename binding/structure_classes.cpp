#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Structure.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "SphericalSurface.hpp"
#include "CuboidalRegion.hpp"
#include "binding_common.hpp"

namespace binding {

void register_structure_classes()
{
    register_structure_class<Structure>("Structure");
    // register_surface_class<Surface>("Surface");
    // register_region_class<Region>("Region");
    register_cuboidal_region_class<CuboidalRegion>("_CuboidalRegion");
    register_planar_surface_class<PlanarSurface>("_PlanarSurface");
    register_spherical_surface_class<SphericalSurface>("_SphericalSurface");
    register_cylindrical_surface_class<CylindricalSurface>("_CylindricalSurface");
}

} // namespace binding

