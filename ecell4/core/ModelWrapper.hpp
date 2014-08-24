#include "Model.hpp"


namespace ecell4
{

class ModelWrapper
{
public:

    typedef Model::species_container_type species_container_type;
    typedef Model::reaction_rule_container_type reaction_rule_container_type;

protected:

    typedef utils::get_mapper_mf<Species::serial_type, Species>::type
        species_attribute_cache_type;
    typedef std::map<Species::serial_type, std::vector<ReactionRule> >
        first_order_reaction_rules_map_type;
    typedef std::map<std::pair<Species::serial_type, Species::serial_type>,
                     std::vector<ReactionRule> >
        second_order_reaction_rules_map_type;

public:

    ModelWrapper(const boost::shared_ptr<Model>& m)
        : model_(m), species_attribute_cache_(), // species_cache_(),
        first_order_reaction_rules_map_(), second_order_reaction_rules_map_()
    {
        ;
    }

    virtual ~ModelWrapper()
    {
        ;
    }

    const boost::shared_ptr<Model>& backend() const
    {
        return model_;
    }

    void initialize()
    {
        // species_cache_.clear();
        species_attribute_cache_.clear();
        first_order_reaction_rules_map_.clear();
        second_order_reaction_rules_map_.clear();
    }

    Integer apply(const Species& pttrn, const Species& sp)
    {
        return model_->apply(pttrn, sp);
    }

    std::vector<ReactionRule> apply(
        const ReactionRule& rr,
        const ReactionRule::reactant_container_type& reactants)
    {
        return model_->apply(rr, reactants);
    }

    Species apply_species_attributes(const Species& sp)
    {
        species_attribute_cache_type::const_iterator
            itr(species_attribute_cache_.find(sp.serial()));
        if (itr != species_attribute_cache_.end())
        {
            return (*itr).second;
        }

        Species retval(model_->apply_species_attributes(sp));
        species_attribute_cache_[sp.serial()] = retval;
        return retval;
    }

    std::vector<ReactionRule> query_reaction_rules(const Species& sp)
    {
        first_order_reaction_rules_map_type::const_iterator
            i(first_order_reaction_rules_map_.find(sp.serial()));
        if (i != first_order_reaction_rules_map_.end())
        {
            return (*i).second;
        }

        std::vector<ReactionRule> retval(model_->query_reaction_rules(sp));
        first_order_reaction_rules_map_.insert(std::make_pair(sp.serial(), retval));
        return retval;
    }

    std::vector<ReactionRule> query_reaction_rules(
        const Species& sp1, const Species& sp2)
    {
        const std::pair<Species::serial_type, Species::serial_type>
            key(sp1.serial() < sp2.serial()?
                std::make_pair(sp1.serial(), sp2.serial()):
                std::make_pair(sp2.serial(), sp1.serial()));
        second_order_reaction_rules_map_type::const_iterator
            i(second_order_reaction_rules_map_.find(key));
        if (i != second_order_reaction_rules_map_.end())
        {
            return (*i).second;
        }

        std::vector<ReactionRule> retval(model_->query_reaction_rules(sp1, sp2));
        second_order_reaction_rules_map_.insert(std::make_pair(key, retval));
        return retval;
    }

    // const std::vector<Species> list_species()
    // {
    //     if (species_cache_.size() == 0)
    //     {
    //         species_cache_ = model_->list_species();
    //     }
    //     return species_cache_;
    // }

protected:

    boost::shared_ptr<Model> model_;

    // species_container_type species_cache_;
    species_attribute_cache_type species_attribute_cache_;
    first_order_reaction_rules_map_type first_order_reaction_rules_map_;
    second_order_reaction_rules_map_type second_order_reaction_rules_map_;
};

} // ecell4
