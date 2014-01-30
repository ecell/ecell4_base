#ifndef __ECELL4_EGFRD_NETWORK_RULES_ADAPTER
#define __ECELL4_EGFRD_NETWORK_RULES_ADAPTER

#include <map>
#include <vector>
#include <boost/scoped_ptr.hpp>

#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>

// cf. epdp/NetworkRulesWrapper.hpp
// NetworkRulesAdapter will substitute for NetworkRulesWrapper 
//  which is instanciated in ParticleSimulatorTraitsBase.
//
//  This class is called via query_reaction_rule function by EGFRDSimulator implemented in epdp, 
//  then, this class translates the query into ecell4::Species or ReactionRule object and 
//  consult ecell4::NetworkModel class.
//template <typename T_, typename Trri_>

template <typename Trri_> 
class NetworkRulesAdapter 
{
public:
    typedef ecell4::NetworkModel backend_type;    // will not be used.
    typedef Trri_ reaction_rule_type;
    typedef typename reaction_rule_type::species_id_type species_id_type;
    typedef std::vector<reaction_rule_type> reaction_rule_vector;
    typedef reaction_rule_vector reaction_rules;
    typedef std::map<species_id_type, reaction_rule_vector> first_order_reaction_rule_vector_map;
    typedef std::map<std::pair<species_id_type, species_id_type>, reaction_rule_vector> second_order_reaction_rule_vector_map;

public:
    reaction_rule_vector const& query_reaction_rule(species_id_type const& r1) const
    {
        // XXX return empty vector.

        // XXX
        //std::vector<ecell4::ReactionRule> ecell4_reaction_rules( this->ecell4_nw_model_->query_reaction_rules(/*XXX*/) );
        return this->empty_reaction_rule_vector_;
    }
    reaction_rule_vector const& query_reaction_rule(
            species_id_type const& r1, species_id_type const& r2) const
    {
        // XXX return empty vector.
        return this->empty_reaction_rule_vector_;
    }

    NetworkRulesAdapter(boost::shared_ptr<ecell4::NetworkModel> ecell4_nw_model):
        ecell4_nw_model_(ecell4_nw_model){;}

private:
    boost::shared_ptr<ecell4::NetworkModel> ecell4_nw_model_;
    mutable first_order_reaction_rule_vector_map first_order_cache_;
    mutable second_order_reaction_rule_vector_map second_order_cache_;
    reaction_rule_vector empty_reaction_rule_vector_;
};



#endif  // __ECELL4_EGFRD_NETWORK_RULES_ADAPTER
