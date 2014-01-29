#ifndef __ECELL4_EGFRD_NETWORK_RULES_ADAPTER
#define __ECELL4_EGFRD_NETWORK_RULES_ADAPTER

#include <map>
#include <vector>
#include <boost/scoped_ptr.hpp>

namespace ecell4 {

namespace egfrd {

// cf. epdp/NetworkRulesWrapper.hpp
// This class will substitute for NetworkRulesWrapper 
//  which is instanciated in ParticleSimulatorTraitsBase.
template <typename T_, typename Trri_>  // will be removed
class NetworkRulesAdapter {
public:
    typedef T_ backend_type;    // will not be used.
    typedef Trri_ reaction_rule_type;
    typedef typename reaction_rule_type::species_id_type species_id_type;
    typedef std::vector<reaction_rule_type> reaction_rule_vector;
    typedef reaction_rule_vector reaction_rules;
    

public:
    reaction_rule_vector const& query_reaction_rule(species_id_type const& r1) const
    {
        reaction_rule_vector retval;
        return retval;
    }
    reaction_rule_vector const& query_reaction_rule(
            species_id_type const& r1, species_id_type const& r2) const
    {
        reaction_rule_vector retval;
        return retval;
    }

    NetworkRulesAdapter(boost::shared_ptr<ecell4::NetworkModel> ecell4_nw_model)
    {;}
};

}   // egfrd

}   // ecell4

#endif  // __ECELL4_EGFRD_NETWORK_RULES_ADAPTER
