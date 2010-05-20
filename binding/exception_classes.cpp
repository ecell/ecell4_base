#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "peer/wrappers/exception/exception_wrapper.hpp"

namespace binding {

void register_exception_classes()
{
    peer::wrappers::exception_wrapper<NotFound, peer::wrappers::py_exc_traits<&PyExc_LookupError> >::__register_class("NotFound");
    peer::wrappers::exception_wrapper<AlreadyExists, peer::wrappers::py_exc_traits<&PyExc_StandardError> >::__register_class("AlreadyExists");
    peer::wrappers::exception_wrapper<IllegalState, peer::wrappers::py_exc_traits<&PyExc_StandardError> >::__register_class("IllegalState");
}

} // namespace binding
