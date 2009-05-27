#ifndef REACTION_RULE_HPP
#define REACTION_RULE_HPP

#include <set>
#include <boost/bind.hpp>

#include "Defs.hpp"

typedef unsigned int reaction_rule_id;

class network_rules;
class species_type;

class reaction_rule
{
    friend class network_rules;
private:
    typedef std::set<species_type*> species_type_set;

public:
    typedef species_type_set::const_iterator species_type_iterator;
    typedef boost::iterator_range<species_type_iterator> species_type_range;

public:
    reaction_rule_id id() const
    {
        return id_;
    }

    void add_reactant(species_type* s)
    {
        reactants_.insert(s);
    }

    species_type_range get_reactants()
    {
        return reactants_;
    }

    void add_product(species_type* s)
    {
        products_.insert(s);
    }

    species_type_range get_products()
    {
        return products_;
    }

    class network_rules* network_rules() const
    {
        return network_rules_;
    }

    Real k() const
    {
        return k_;
    }

    Real& k()
    {
        return k_;
    }

    reaction_rule() {}

protected:
    class network_rules*& network_rules()
    {
        return network_rules_;
    }

    reaction_rule_id& id()
    {
        return id_;
    }

private:
    class network_rules* network_rules_;
    species_type_set reactants_;
    species_type_set products_;
    reaction_rule_id id_;
    Real k_;
};

template<typename T1_, typename T2_>
reaction_rule* new_reaction_rule(T1_ const& reactants, T2_ const& products, Real k)
{
    reaction_rule* retval = new reaction_rule();
    retval->k() = k;
    std::for_each(boost::begin(reactants), boost::end(reactants),
            boost::bind(&reaction_rule::add_reactant, retval, _1));
    std::for_each(boost::begin(products), boost::end(products),
            boost::bind(&reaction_rule::add_product, retval, _1));
    return retval;
}

#endif /* REACTION_RULE_HPP */
