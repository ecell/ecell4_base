#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include "peer/util/to_native_converter.hpp"
#include "binding_common.hpp"
#include "../ConsoleAppender.hpp"

namespace binding {

boost::python::objects::class_base
register_console_appender_class(char const* name)
{
    using namespace boost::python;
    typedef ConsoleAppender impl_type;

    return class_<impl_type, bases<impl_type::base_type>,
                  boost::shared_ptr<impl_type>, boost::noncopyable>(name, init<>())
        ;
}

} // namespace binding
