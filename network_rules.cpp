#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "network_rules.hpp"
#include "reaction_rule.hpp"
#include "species_type.hpp"

network_rules::network_rules(): last_id_(0)
{
}

void network_rules::add_reaction_rule(reaction_rule* r)
{
    r->network_rules() = this;
    r->id() = last_id_++;
    reaction_rules_.insert(std::make_pair(r->id(), r));
}

