#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "PythonLogger.hpp"
#include "binding_common.hpp"

namespace binding {

void register_python_logger_classes()
{
    register_logger_factory_class("_LoggerFactory");
    register_python_logger_factory_class("PythonLoggerFactory");
}

} // namespace binding
