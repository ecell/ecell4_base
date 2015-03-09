#ifndef __ECELL4_ODE_ODE_NETWORK_MODEL_HPP
#define __ECELL4_ODE_ODE_NETWORK_MODEL_HPP

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
    ODENetworkModel(const boost::shared_ptr<ecell4::NetworkModel> model);
    ~ODENetworkModel();

    void update_model();
    bool has_model() const
    {
        return !(this->expanded_.expired());
    }
    boost::shared_ptr<NetworkModel> get_networkmodel() const
    {
        return this->expanded_.lock();
    }

    const ode_reaction_rule_container_type& ode_reaction_rules() const
    {
        return ode_reaction_rules_;
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
private:
    bool convert_from_networkmodel(const boost::shared_ptr<ecell4::NetworkModel> model);

protected:
    species_container_type species_attributes_;
    ode_reaction_rule_container_type ode_reaction_rules_;

    first_order_reaction_rules_map_type first_order_reaction_rules_map_;
    second_order_reaction_rules_map_type second_order_reaction_rules_map_;

    boost::weak_ptr<NetworkModel> expanded_;
};

}   // ode

}   // ecell4

#endif
