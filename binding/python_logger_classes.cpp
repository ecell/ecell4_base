#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "PythonAppender.hpp"
#include "Logger.hpp"
#include "LoggerManager.hpp"
#include "LogAppender.hpp"
#include "ConsoleAppender.hpp"

namespace binding {

void register_python_logger_classes()
{
    register_logger_handler_class("CppLoggerHandler");
    register_logger_class("Logger");
    register_logger_manager_class("LoggerManager");
    register_log_appender_class("LogAppender");
    register_console_appender_class("ConsoleAppender");
    register_python_appender_class("PythonAppender");
    register_logger_level_enum("LogLevel");
}

} // namespace binding
