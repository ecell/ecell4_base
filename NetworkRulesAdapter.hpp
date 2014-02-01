#ifndef __ECELL4_EGFRD_NETWORK_RULES_ADAPTER
#define __ECELL4_EGFRD_NETWORK_RULES_ADAPTER

#include <map>
#include <numeric>
#include <vector>
#include <boost/scoped_ptr.hpp>
#include "twofold_container.hpp"
#include "ReactionRuleInfo.hpp"
#include "exceptions.hpp"

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
        typename first_order_reaction_rule_vector_map::const_iterator i(
                this->first_order_cache_.find(r1));
        if (i == this->first_order_cache_.end())
        {
            ecell4::NetworkModel::reaction_rule_container_type 
                reaction_rules_at_ecell4( this->ecell4_nw_model_->query_reaction_rules(ecell4::Species(r1)) );

            std::pair<typename first_order_reaction_rule_vector_map::iterator, bool> 
                x(this->first_order_cache_.insert( std::make_pair(r1, reaction_rule_vector()) ));
                
            for (std::vector<ecell4::ReactionRule>::iterator it( reaction_rules_at_ecell4.begin() );
                    it != reaction_rules_at_ecell4.end();
                    it++ ) 
            {
                //this->temporary_reaction_rule_vector_.push_back( convert_reaction_rule_type(*it) );
                x.first->second.push_back(convert_reaction_rule_type( *it ));
            }
            return x.first->second;
        }
        return i->second;
    }

    reaction_rule_vector const& query_reaction_rule(
            species_id_type const& r1, species_id_type const& r2) const
    {
        typename second_order_reaction_rule_vector_map::const_iterator i(
                this->second_order_cache_.find(std::make_pair(r1, r2)) );
        //asm volatile ("int3");
        if (i == this->second_order_cache_.end())
        {
            ecell4::NetworkModel::reaction_rule_container_type 
                reaction_rules_at_ecell4( this->ecell4_nw_model_->query_reaction_rules(ecell4::Species(r1)) );

            std::pair<typename second_order_reaction_rule_vector_map::iterator, bool>
                x(this->second_order_cache_.insert( std::make_pair(std::make_pair(r1, r2), reaction_rule_vector()) ));

            for (std::vector<ecell4::ReactionRule>::iterator it( reaction_rules_at_ecell4.begin() );
                    it != reaction_rules_at_ecell4.end();
                    it++ ) 
            {
                x.first->second.push_back(convert_reaction_rule_type( *it ));
            }
            return x.first->second;
        }
        return i->second;
    }

    NetworkRulesAdapter(boost::shared_ptr<ecell4::NetworkModel> ecell4_nw_model):
        ecell4_nw_model_(ecell4_nw_model) 
    {;}

protected:
    inline reaction_rule_type convert_reaction_rule_type(const ecell4::ReactionRule& rr) const
    {
        //typedef twofold_container<species_id_type> reactants_container_type;
        typedef typename reaction_rule_type::rate_type rate_type;

        reaction_rule_type retval;
        std::vector<species_id_type> products;
        std::vector<species_id_type> reactants;
        rate_type rate;
        try {
            rate = boost::lexical_cast<rate_type>(rr.k());
        } 
        catch (boost::bad_lexical_cast &) 
        {
            if (rr.k() == double(HUGE_VAL)) {
                rate = std::numeric_limits<rate_type>::infinity();
            }
            else {
                throw;
            }
        }
        for (ecell4::ReactionRule::product_container_type::const_iterator
                 j(rr.products().begin()); j != rr.products().end(); ++j)
        {
            products.push_back( j->name() );
        }

        ecell4::ReactionRule::reactant_container_type::const_iterator
            r(rr.reactants().begin());
        switch (rr.reactants().size())
        {
        case 1:
            {
                //const ::SpeciesTypeID sid1(find(*r));
                const species_id_type sid1(r->name());
                reactants.push_back(sid1);
                return reaction_rule_type(rr.id(), rate, reactants, products);
            }
            break;
        case 2:
            {
                //const ::SpeciesTypeID sid1(find(*r));
                const species_id_type sid1(r->name());
                reactants.push_back(sid1);
                ++r;
                //const ::SpeciesTypeID sid2(find(*r));
                const species_id_type sid2(r->name());
                reactants.push_back(sid2);
                return reaction_rule_type(rr.id(), rate, reactants, products);
            }
            break;
        default:
            throw illegal_state("the number of reactants must be 1 or 2.");
            break;
        }
        return reaction_rule_type(); // never get here
    }


private:
    boost::shared_ptr<ecell4::NetworkModel> ecell4_nw_model_;
    mutable first_order_reaction_rule_vector_map first_order_cache_;
    mutable second_order_reaction_rule_vector_map second_order_cache_;
};



#endif  // __ECELL4_EGFRD_NETWORK_RULES_ADAPTER
