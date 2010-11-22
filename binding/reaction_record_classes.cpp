#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include "ReactionRecord.hpp"
#include "reaction_recorder_converter.hpp"
#include "binding_common.hpp"

namespace binding {

void register_reaction_record_classes()
{
    register_reaction_record_class<ReactionRecord>("ReactionRecord");
    register_reaction_recorder_converter<ReactionRecorder>();
}

} // namespace binding
