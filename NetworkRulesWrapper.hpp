#ifndef NETWORK_RULES_WRAPPER_HPP
#define NETWORK_RULES_WRAPPER_HPP

#include <map>
#include <vector>
#include <boost/scoped_ptr.hpp>
#include <boost/bind.hpp>
#include "twofold_container.hpp"
#include "utils/range.hpp"
#include "ReactionRuleInfo.hpp"

template<typename T_, typename Trri_>
class NetworkRulesWrapper
{
public:
    typedef T_ backend_type;
    typedef Trri_ reaction_rule_type;
    typedef typename reaction_rule_type::species_id_type species_id_type;
    typedef std::vector<reaction_rule_type> reaction_rule_vector;
    typedef reaction_rule_vector reaction_rules;
    typedef std::map<species_id_type, reaction_rule_vector> first_order_reaction_rule_vector_map;
    typedef std::map<std::pair<species_id_type, species_id_type>, reaction_rule_vector> second_order_reaction_rule_vector_map;

public:
    reaction_rule_vector const& query_reaction_rule(species_id_type const& r1) const
    {
        typename first_order_reaction_rule_vector_map::const_iterator i(
            first_order_cache_.find(r1));
        if (i == first_order_cache_.end())
        {
            std::pair<
                typename first_order_reaction_rule_vector_map::iterator,
                bool> x(first_order_cache_.insert(
                    std::make_pair(r1, reaction_rule_vector())));
            boost::scoped_ptr<typename backend_type::reaction_rule_generator>
                gen(backend_.query_reaction_rule(r1));
            if (gen)
            {
                while (::valid(*gen))
                {
                    typename backend_type::reaction_rule_type const r((*gen)());
                    (*x.first).second.push_back(reaction_rule_type(
                        r.id(),
                        boost::lexical_cast<
                            typename reaction_rule_type::rate_type>(r["k"]),
                        r.get_reactants(),
                        r.get_products()));
                }
            }
            return (*x.first).second;
        }
        return (*i).second;
    }

    reaction_rule_vector const& query_reaction_rule(
            species_id_type const& r1, species_id_type const& r2) const
    {
        typename second_order_reaction_rule_vector_map::const_iterator i(
            second_order_cache_.find(std::make_pair(r1, r2)));
        if (i == second_order_cache_.end())
        {
            std::pair<
                typename second_order_reaction_rule_vector_map::iterator,
                bool> x(second_order_cache_.insert(
                    std::make_pair(std::make_pair(r1, r2),
                                   reaction_rule_vector())));
            boost::scoped_ptr<typename backend_type::reaction_rule_generator>
                gen(backend_.query_reaction_rule(r1, r2));
            if (gen)
            {
                while (::valid(*gen))
                {
                    typename backend_type::reaction_rule_type const r((*gen)());
                    (*x.first).second.push_back(reaction_rule_type(
                        r.id(),
                        boost::lexical_cast<
                            typename reaction_rule_type::rate_type>(r["k"]),
                        r.get_reactants(),
                        r.get_products()));
                }
            }
            return (*x.first).second;
        }
        return (*i).second;
    }

    NetworkRulesWrapper(backend_type const& backend): backend_(backend) {}

private:
    mutable first_order_reaction_rule_vector_map first_order_cache_;
    mutable second_order_reaction_rule_vector_map second_order_cache_;
    backend_type const& backend_;
};

#endif /* NETWORK_RULES_WRAPPER_HPP */
