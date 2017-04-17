#ifndef ECELL4_ODE_ODE_NETWORK_MODEL_HPP
#define ECELL4_ODE_ODE_NETWORK_MODEL_HPP

#include <map>
#include <set>
#include <algorithm>
#include <iterator>

// XXX debug;
#include <iostream>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "ODEReactionRule.hpp"


namespace ecell4
{
namespace ode
{

class ODENetworkModel
    //: public ecell4::Model
{
public:
    // derive from Core class' typedefs
    typedef Model abstract_model_type;
    typedef NetworkModel network_model_type;
    typedef abstract_model_type::species_container_type species_container_type;
    
    // this class specific typedefs.
    typedef ODEReactionRule ode_reaction_rule_type;
    typedef std::vector<ODEReactionRule> ode_reaction_rule_container_type;
protected:
    typedef std::map<Species::serial_type, 
                     std::vector<ode_reaction_rule_container_type::size_type> >
        first_order_reaction_rules_map_type;
    typedef std::map<std::pair<Species::serial_type, Species::serial_type>, 
                     std::vector<ode_reaction_rule_container_type::size_type> >
        second_order_reaction_rules_map_type;
public:
    ODENetworkModel();
    ODENetworkModel(const boost::shared_ptr<ecell4::Model> model);
    ~ODENetworkModel();

    void update_model();
    bool has_network_model() const
    {
        return !(this->expanded_.expired());
    }

    const std::vector<Species> list_species() const
    {
        std::vector<Species> retval;
        const ode_reaction_rule_container_type &rrs(ode_reaction_rules());
        for(ode_reaction_rule_container_type::const_iterator i(rrs.begin());
            i != rrs.end(); i++)
        {
            const ODEReactionRule::reactant_container_type &reactants((*i).reactants());
            const ODEReactionRule::product_container_type &products((*i).products());
            std::copy(reactants.begin(), reactants.end(), std::back_inserter(retval));
            std::copy(products.begin(), products.end(), std::back_inserter(retval));
        }
        std::sort(retval.begin(), retval.end());
        retval.erase(std::unique(retval.begin(), retval.end()), retval.end());
        return retval;
    }
    boost::shared_ptr<Model> get_networkmodel() const
    {
        return this->expanded_.lock();
    }

    const ode_reaction_rule_container_type& ode_reaction_rules() const
    {
        return ode_reaction_rules_;
    }
    inline const ode_reaction_rule_container_type& reaction_rules() const
    {
        return ode_reaction_rules();
    }
    const species_container_type& species_attributes() const
    {
        return species_attributes_;
    }
    const Integer num_reaction_rules() const
    {
        return ode_reaction_rules_.size();
    }

    void dump_reactions() const
    {
        for(ode_reaction_rule_container_type::const_iterator it(ode_reaction_rules_.begin());
                it != ode_reaction_rules_.end(); it++)
        {
            std::cout << it->as_string() << std::endl;
        }
    }
    void add_reaction_rule(ODEReactionRule ode_rr)
    {
        this->ode_reaction_rules_.push_back(ode_rr);
    }
    void add_reaction_rule(ReactionRule ode_rr)
    {
        this->ode_reaction_rules_.push_back(ODEReactionRule(ode_rr));
    }

    void add_reaction_rules(const std::vector<ODEReactionRule>& rrs)
    {
        for (std::vector<ODEReactionRule>::const_iterator i(rrs.begin());
            i != rrs.end(); ++i)
        {
            add_reaction_rule(*i);
        }
    }

    void add_reaction_rules(const std::vector<ReactionRule>& rrs)
    {
        for (std::vector<ReactionRule>::const_iterator i(rrs.begin());
            i != rrs.end(); ++i)
        {
            add_reaction_rule(*i);
        }
    }

    void add_species_attribute(const Species &sp)
    {
        if (has_species_attribute(sp))
        {
            throw AlreadyExists("species already exista");
        }
        species_attributes_.push_back(sp);
    }

    bool has_species_attribute(const Species& sp) const
    {
        species_container_type::const_iterator i(
            std::find(species_attributes_.begin(), species_attributes_.end(), sp));
        return (i != species_attributes_.end());
    }
private:
    bool convert_from_networkmodel(const boost::shared_ptr<ecell4::Model> model);

protected:
    species_container_type species_attributes_;
    ode_reaction_rule_container_type ode_reaction_rules_;

    first_order_reaction_rules_map_type first_order_reaction_rules_map_;
    second_order_reaction_rules_map_type second_order_reaction_rules_map_;

    boost::weak_ptr<Model> expanded_;
};

}   // ode

}   // ecell4

#endif
