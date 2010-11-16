#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "PythonLogger.hpp"
#include "binding_common.hpp"

namespace binding {

void register_python_logger_classes()
{
    register_logger_handler_class("CppLoggerHandler");
    register_logger_factory_class("LoggerFactory");
    register_null_logger_factory_class("NullLoggerFactory");
    register_python_logger_factory_class("PythonLoggerFactory");
    register_logger_level_enum("LogLevel");
}

} // namespace binding
