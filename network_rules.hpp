#ifndef NETWORK_RULES_HPP
#define NETWORK_RULES_HPP

#include <map>
#include <boost/range/iterator_range.hpp>
#include "reaction_rule.hpp"

class network_rules
{
    typedef std::map<reaction_rule_id, reaction_rule*> reaction_rule_map;

public:
    typedef reaction_rule_map::const_iterator reaction_rule_iterator;
    typedef boost::iterator_range<reaction_rule_iterator> reaction_rule_range;

public:
    void add_reaction_rule(reaction_rule* reaction_rule);

    network_rules();

private:
    reaction_rule_map reaction_rules_;
    reaction_rule_id last_id_;
};


#endif /* NETWORK_RULES_HPP */
