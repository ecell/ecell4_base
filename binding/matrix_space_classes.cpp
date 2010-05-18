#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "MatrixSpace.hpp"
#include "binding_common.hpp"

namespace binding {

void register_spherical_shell_container_class()
{
    register_matrix_space_class<SphericalShellContainer>(
            "SphericalShellContainer");
}

void register_cylindrical_shell_container_class()
{
    register_matrix_space_class<CylindricalShellContainer>(
            "CylindricalShellContainer");
}

} // namespace binding
