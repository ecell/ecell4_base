#ifndef REACTION_RULE_HPP
#define REACTION_RULE_HPP

#include <vector>
#include <ostream>
#include <algorithm>

#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

#include "Defs.hpp"
#include "utils/get_mapper_mf.hpp"
#include "utils/range_support.hpp"
#include "exceptions.hpp"
#include "utils.hpp"
#include "SpeciesTypeID.hpp"
#include "twofold_container.hpp"
#include "utils/memberwise_compare.hpp"

class NetworkRules;

class ReactionRule
{
public:
    typedef SpeciesTypeID species_type_id_type;

private:
    typedef std::vector<species_type_id_type> species_type_id_vector;
    typedef get_mapper_mf<std::string, std::string>::type string_map_type;

public:
    typedef int identifier_type; 
    typedef string_map_type::const_iterator string_map_iterator;
    typedef boost::iterator_range<string_map_iterator> attributes_range;

    typedef species_type_id_vector::const_iterator species_type_id_iterator;
    typedef boost::iterator_range<species_type_id_iterator> species_type_id_range;
    typedef twofold_container<SpeciesTypeID> Reactants;

public:
    Reactants const& get_reactants() const
    {
        return reactants_;
    }

    void add_product(SpeciesTypeID const& s)
    {
        products_.insert(
            std::lower_bound(products_.begin(), products_.end(), s),
            s);
    }

    species_type_id_range get_products() const
    {
        return products_;
    }

    identifier_type const& id() const
    {
        return id_;
    }

    // package-private
    identifier_type const& set_id(identifier_type const& val) const
    {
        id_ = val;
        return id_;
    }

    std::string const& operator[](std::string const& name) const
    {
        string_map_type::const_iterator i(attrs_.find(name));
        if (i == attrs_.end())
            throw not_found((boost::format("key %s not found") % name).str());
        return (*i).second;
    }

    std::string& operator[](std::string const& name)
    {
        return attrs_[name];
    }

    attributes_range attributes() const
    {
        return attributes_range(attrs_.begin(), attrs_.end());
    }

    ReactionRule()
        : id_(), reactants_() {}

    explicit ReactionRule(Reactants const& _reactants)
        : id_(), reactants_(_reactants) {}

    template<typename Trange_>
    ReactionRule(Reactants const& _reactants, Trange_ const& products)
        : id_(), reactants_(_reactants)
    {
        std::for_each(boost::begin(products), boost::end(products),
                boost::bind(&ReactionRule::add_product, this, _1));
        std::stable_sort(products_.begin(), products_.end());
    }

private:
    mutable identifier_type id_;
    Reactants reactants_;
    species_type_id_vector products_;
    Real k_;
    string_map_type attrs_;
};

inline bool operator<(ReactionRule const& lhs, ReactionRule const& rhs)
{
    int tmp = memberwise_compare(lhs.get_reactants(), rhs.get_reactants());
    if (tmp > 0)
    {
        return false;
    }
    else if (tmp < 0)
    {
        return true;
    }
    return memberwise_compare(lhs.get_products(), rhs.get_products()) < 0;
}

inline bool operator==(ReactionRule const& lhs, ReactionRule const& rhs)
{
    return lhs.get_reactants() == rhs.get_reactants() &&
            memberwise_compare(lhs.get_products(), rhs.get_products()) == 0;
}

inline bool operator!=(ReactionRule const& lhs, ReactionRule const& rhs)
{
    return !(lhs == rhs);
}

template<typename T2_>
inline ReactionRule new_reaction_rule(SpeciesTypeID const& r1, T2_ const& products, Real k)
{
    ReactionRule retval((ReactionRule::Reactants(r1)));
    retval["k"] = boost::lexical_cast<std::string>(k);
    std::for_each(boost::begin(products), boost::end(products),
            boost::bind(&ReactionRule::add_product, &retval, _1));
    return retval;
}

template<typename T2_>
inline ReactionRule new_reaction_rule(SpeciesTypeID const& r1, SpeciesTypeID const& r2, T2_ const& products, Real k)
{
    ReactionRule retval(ReactionRule::Reactants(r1, r2));
    retval["k"] = boost::lexical_cast<std::string>(k);
    std::for_each(boost::begin(products), boost::end(products),
            boost::bind(&ReactionRule::add_product, &retval, _1));
    return retval;
}

inline bool valid(ReactionRule const& r)
{
    return r.get_reactants().size() != 0;
}

template<typename Tchar_, typename Ttraits_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& out, ReactionRule const& r)
{
    bool first;
    out << "ReactionRule(id=" << r.id() << ", reactants={";
    first = true;
    ReactionRule::Reactants const& reactants(r.get_reactants());
    ReactionRule::species_type_id_range products(r.get_products());
    for (typename boost::range_const_iterator<ReactionRule::Reactants>::type
            i(boost::begin(reactants)), e(boost::end(reactants));
         i != e; ++i)
    {
        if (!first)
        {
            out << ", ";
        }
        out << *i;
        first = false;
    }
    out << "}, products={";
    first = true;
    for (typename boost::range_const_iterator<
            ReactionRule::species_type_id_range>::type
                i(boost::begin(products)), e(boost::end(products));
         i != e; ++i)
    {
        if (!first)
        {
            out << ", ";
        }
        out << *i;
        first = false;
    }
    out << "})";
    return out;
}

#endif /* REACTION_RULE_HPP */
