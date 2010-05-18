#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Shell.hpp"
#include "binding_common.hpp"

namespace binding {

void register_spherical_shell_class()
{
    ShellWrapper<SphericalShell>::__register_class("SphericalShell");
}

void register_cylindrical_shell_class()
{
    ShellWrapper<CylindricalShell>::__register_class("CylindricalShell");
}

} // namespace binding
