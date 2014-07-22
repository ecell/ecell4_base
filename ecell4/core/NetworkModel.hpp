#ifndef __ECELL4_NETWORK_MODEL_HPP
#define __ECELL4_NETWORK_MODEL_HPP

// #include "get_mapper_mf.hpp"

#include <map>
#include <set>
#include <algorithm>
#include <iterator>

#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "Model.hpp"

#include "Context.hpp"


namespace ecell4
{

class NetworkModel
    : public Model
{
public:

    typedef std::vector<Species> species_container_type;
    typedef std::vector<ReactionRule> reaction_rule_container_type;

protected:


    typedef std::map<Species::serial_type,
                     std::vector<reaction_rule_container_type::size_type> >
        first_order_reaction_rules_map_type;
    typedef std::map<std::pair<Species::serial_type, Species::serial_type>,
                     std::vector<reaction_rule_container_type::size_type> >
        second_order_reaction_rules_map_type;

    typedef utils::get_mapper_mf<
        Species::serial_type, species_container_type::size_type>::type
        species_attribute_cache_type;

public:

    NetworkModel()
        : dirty_(false), species_attributes_(), species_cache_(), reaction_rules_(),
        first_order_reaction_rules_map_(), second_order_reaction_rules_map_(),
        species_attribute_cache_()
    {
        ;
    }

    virtual ~NetworkModel()
    {
        ;
    }

    // ModelTraits

    std::vector<ReactionRule> query_reaction_rules(const Species& sp) const;
    std::vector<ReactionRule> query_reaction_rules(
        const Species& sp1, const Species& sp2) const;

    Species apply_species_attributes(const Species& sp) // const
    {
        species_attribute_cache_type::const_iterator
            itr(species_attribute_cache_.find(sp.serial()));
        if (itr != species_attribute_cache_.end())
        {
            Species retval(sp);
            retval.set_attributes(species_attributes_[(*itr).second]);
            return retval;
        }

        for (species_container_type::iterator
            i(species_attributes_.begin()); i != species_attributes_.end(); ++i)
        {
            // if (sp == *i)
            if (spmatch(*i, sp))
            {
                Species retval(sp);
                retval.set_attributes(*i);
                species_attribute_cache_[sp.serial()] =
                    std::distance(species_attributes_.begin(), i);
                return retval;
            }
        }
        return sp;
    }

    // NetworkModelTraits

    void add_species_attribute(const Species& sp);
    bool has_species_attribute(const Species& sp) const;
    void remove_species_attribute(const Species& sp);

    void add_reaction_rule(const ReactionRule& rr);
    void remove_reaction_rule(const ReactionRule& rr);
    bool has_reaction_rule(const ReactionRule& rr) const;

    // Optional functions

    const std::vector<Species> list_species()
    {
        initialize();
        return species_cache_;
    }

    Species create_species(const std::string& name) // const
    {
        return apply_species_attributes(Species(name));
    }

    const reaction_rule_container_type& reaction_rules() const
    {
        return reaction_rules_;
    }

    const species_container_type& species_attributes() const
    {
        return species_attributes_;
    }

    const Integer num_reaction_rules() const
    {
        return reaction_rules_.size();
    }

protected:

    void initialize();

protected:

    bool dirty_;
    species_container_type species_attributes_, species_cache_;
    reaction_rule_container_type reaction_rules_;

    first_order_reaction_rules_map_type first_order_reaction_rules_map_;
    second_order_reaction_rules_map_type second_order_reaction_rules_map_;

    species_attribute_cache_type species_attribute_cache_;
};

} // ecell4

#endif /* __ECELL4_NETWORK_MODEL_HPP */
