#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Model.hpp"
#include "binding_common.hpp"

namespace binding {

void register_model_class()
{
    register_model_class<Model>("Model");
}

} // namesapce binding
