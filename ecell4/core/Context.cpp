#include "Context.hpp"
#include <string>
#include <sstream>


namespace ecell4
{

namespace context
{

std::string itos(unsigned int val)
{
    std::stringstream ss;
    ss << val;
    return ss.str();
}

unsigned int concatenate_units(std::vector<UnitSpecies>& units1, const Species& sp, const unsigned int bond_stride)
{
    const std::vector<UnitSpecies> units2 = sp.units();
    units1.reserve(units1.size() + units2.size());

    std::unordered_map<std::string, std::string> bond_cache;

    for (Species::container_type::const_iterator j(units2.begin());
        j != units2.end(); ++j)
    {
        units1.push_back(*j);

        for (UnitSpecies::container_type::const_iterator
            k((*j).begin()); k != (*j).end(); ++k)
        {
            const std::string& bond((*k).second.second);
            if (bond != "" && !is_wildcard(bond))
            {
                std::unordered_map<std::string, std::string>::const_iterator
                    it = bond_cache.find(bond);
                if (it == bond_cache.end())
                {
                    const std::string newbond = itos(bond_stride + 1 + bond_cache.size());
                    bond_cache.insert(std::make_pair(bond, newbond));
                    units1.back().at(std::distance((*j).begin(), k)).second.second = newbond;
                }
                else
                {
                    units1.back().at(std::distance((*j).begin(), k)).second.second = (*it).second;
                }
            }
        }
    }
    return bond_cache.size();
}

unsigned int tag_units(
    std::vector<unsigned int>& groups,
    const unsigned int& group_id,
    const unsigned int& idx,
    const std::vector<std::vector<std::vector<UnitSpecies>::size_type> >& connections,
    const unsigned int& notyet)
{
    if (groups[idx] != notyet)
    {
        // assert(groups[idx] < group_id);
        return group_id;
    }

    groups[idx] = group_id;

    for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
        i(connections[idx].begin()); i != connections[idx].end(); ++i)
    {
        tag_units(groups, group_id, *i, connections, notyet);
        // const unsigned int _gid = tag_units(groups, group_id, *i, connections, notyet);
        // assert(_gid == group_id);
    }

    return group_id + 1;
}

bool is_correspondent(const UnitSpecies& usp1, const UnitSpecies& usp2)
{
    if (usp1.name() != usp2.name() || usp1.num_sites() != usp2.num_sites())
    {
        return false;
    }

    for (UnitSpecies::container_type::const_iterator
        i(usp1.begin()), j(usp2.begin()); i != usp1.end() && j != usp2.end(); ++i, ++j)
    {
        if ((*i).first != (*j).first)
        {
            return false;
        }
    }
    return true;
}

// struct species_formatter
// {
//     typedef std::vector<UnitSpecies>::size_type size_type;
// 
//     struct bond_type
//     {
//         std::string name;
//         size_type one;
//         size_type another;
// 
//         bond_type()
//         {}
// 
//         bond_type(const std::string& name, const size_type one, const size_type another)
//             : name(name), one(one), another(another)
//         {}
//     };
// 
//     typedef std::unordered_map<std::string, bond_type> bond_container_type;
// 
//     species_formatter(const std::vector<UnitSpecies>& units)
//         : units(units)
//     {
//         bonds.clear();
//         for (size_type idx = 0; idx != units.size(); ++idx)
//         {
//             const UnitSpecies& unit(units.at(idx));
//             for (UnitSpecies::container_type::const_iterator i(unit.begin());
//                  i != unit.end(); ++i)
//             {
//                 const std::string& bond = (*i).second.second;
// 
//                 if (is_empty(bond) || is_wildcard(bond))
//                 {
//                     continue;
//                 }
// 
//                 bond_container_type::iterator it = bonds.find(bond);
//                 if (it == bonds.end())
//                 {
//                     bonds.insert(std::make_pair(bond, bond_type(bond, idx, units.size())));
//                 }
//                 else
//                 {
//                     (*it).second.another = idx;
//                 }
//             }
//         }
//     }
// 
//     size_type get_partner(const std::string bond, const size_type one) const
//     {
//         bond_container_type::const_iterator i = bonds.find(bond);
//         assert(i != bonds.end());
//         const bond_type& b = (*i).second;
//         if (b.one == one)
//         {
//             assert(b.another != units.size());
//             return b.another;
//         }
//         else
//         {
//             assert(b.one != units.size());
//             return b.one;
//         }
//     }
// 
//     // std::vector<size_type> sort(const size_type start) const
//     // {
//     //     std::vector<size_type> permutation(units.size(), units.size());
// 
//     //     size_type stride = sort(permutation, start, 0);
//     //     for (size_type idx = 0; idx != units.size(); ++idx)
//     //     {
//     //         stride = sort(permutation, idx, stride);
//     //     }
// 
//     //     assert(stride == units.size());
//     //     return permutation;
//     // }
// 
//     size_type sort(std::vector<size_type>& permutation, const size_type pos, const size_type idx) const
//     {
//         if (permutation[pos] != units.size())
//         {
//             return idx;
//         }
// 
//         permutation[pos] = idx;
//         size_type stride = idx + 1;
//         const UnitSpecies& unit = units.at(pos);
//         for (UnitSpecies::container_type::const_iterator i(unit.begin());
//              i != unit.end(); ++i)
//         {
//             const std::string& bond = (*i).second.second;
// 
//             if (is_empty(bond) || is_wildcard(bond))
//             {
//                 continue;
//             }
// 
//             stride = sort(permutation, get_partner(bond, pos), stride);
//         }
//         return stride;
//     }
// 
//     Species collect(const std::vector<size_type>& permutation) const
//     {
//         typedef std::unordered_map<std::string, std::string> mapper_type;
//         mapper_type bond_names;
// 
//         std::vector<UnitSpecies> res;
//         for (size_type i = 0; i != permutation.size(); ++i)
//         {
//             const size_type& j = permutation.at(i);
//             if (j == units.size())
//             {
//                 continue;
//             }
//             else if (j >= res.size())
//             {
//                 res.resize(j + 1);
//             }
//             res[j] = units[i];
//         }
// 
//         unsigned int stride = 0;
//         Species sp;
//         for (std::vector<UnitSpecies>::iterator i = res.begin();
//             i != res.end(); ++i)
//         {
//             UnitSpecies& unit = (*i);
// 
//             for (UnitSpecies::container_type::iterator k(unit.begin());
//                  k != unit.end(); ++k)
//             {
//                 std::string& bond = (*k).second.second;
// 
//                 if (is_empty(bond) || is_wildcard(bond))
//                 {
//                     continue;
//                 }
// 
//                 mapper_type::const_iterator
//                     it = bond_names.find(bond);
//                 if (it == bond_names.end())
//                 {
//                     const std::string new_bond_name = itos(++stride);
//                     bond_names.insert(std::make_pair(bond, new_bond_name));
//                     bond = new_bond_name;
//                 }
//                 else
//                 {
//                     bond = (*it).second;
//                 }
//             }
// 
//             sp.add_unit(unit);
//         }
// 
//         return sp;
//     }
// 
//     const std::vector<UnitSpecies>& units;
//     bond_container_type bonds;
// };

class species_structure
{
public:

    // typedef Species::container_type::size_type index_type;
    typedef unsigned int index_type;
    typedef std::pair<index_type, std::string> site_type;
    typedef std::unordered_map<std::string, std::vector<site_type>>
        connection_container_type;

public:

    species_structure(const Species& sp)
        : root_(sp.units())
    {
        initialize();
    }

    const std::vector<UnitSpecies>& units() const
    {
        return root_;
    }

    void initialize()
    {
        connections_.clear();
        for (index_type idx(0); idx < root_.size(); ++idx)
        {
            const UnitSpecies usp(root_.at(idx));
            for (UnitSpecies::container_type::const_iterator i(usp.begin());
                 i != usp.end(); ++i)
            {
                if ((*i).second.second == "" || is_wildcard((*i).second.second))
                {
                    continue;
                }

                if (connections_.find((*i).second.second) == connections_.end())
                {
                    connections_.insert(std::make_pair(
                        (*i).second.second, std::vector<site_type>()));
                }
                connections_[(*i).second.second].push_back(
                    std::make_pair(idx, (*i).first));
            }
        }
    }

    int compare(const index_type& val1, const index_type& val2)
    {
        if (val1 == val2)
        {
            return 0;
        }

        const std::pair<index_type, index_type> pair_key((val1 < val2)?
            std::make_pair(val1, val2) : std::make_pair(val1, val2));
        if (std::binary_search(ignores_.begin(), ignores_.end(), pair_key))
        {
            return 0;
        }

        const UnitSpecies& lhs(root_.at(val1));
        const UnitSpecies& rhs(root_.at(val2));

        if (lhs.name() != rhs.name())
        {
            return (lhs.name() < rhs.name()? 1 : -1);
        }

        UnitSpecies::container_type::const_iterator
            i(lhs.begin()), j(rhs.begin());
        while (i != lhs.end() && j != rhs.end())
        {
            if ((*i).first != (*j).first)
            {
                // std::cout << "[1] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ") -> " << (*i).first
                //     << " < " << (*j).first << std::endl;
                return ((*i).first < (*j).first? 1 : -1);
            }
            else if ((*i).second.first != (*j).second.first)
            {
                // std::cout << "[2] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ")" << std::endl;
                return ((*i).second.first < (*j).second.first? 1 : -1);
            }
            else if (((*i).second.second == "") != ((*j).second.second == ""))
            {
                // std::cout << "[3] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ") -> '"
                //     << (*i).second.second << "' < '" << (*j).second.second
                //     << "'" << std::endl;
                return ((*i).second.second == ""? 1 : -1);
            }

            ++i;
            ++j;
        }

        if (lhs.num_sites() != rhs.num_sites())
        {
            return (lhs.num_sites() < rhs.num_sites()? 1 : -1);
        }

        ignores_.insert(
            std::lower_bound(ignores_.begin(), ignores_.end(), pair_key),
            pair_key);
        i = lhs.begin();
        j = rhs.begin();
        while (i != lhs.end() && j != rhs.end())
        {
            if ((*i).second.second != "" && (*i).second.second != "")
            {
                const std::vector<site_type>&
                    pair1(connections_[(*i).second.second]);
                const std::vector<site_type>&
                    pair2(connections_[(*j).second.second]);
                const site_type& target1(
                    (pair1[0].first == val1 && pair1[0].second == (*i).first)?
                    pair1[1] : pair1[0]);
                const site_type& target2(
                    (pair2[0].first == val2 && pair2[0].second == (*j).first)?
                    pair2[1] : pair2[0]);
                if (target1.second != target2.second)
                {
                    ignores_.pop_back();
                    return (target1.second < target2.second? 1 : -1);
                }

                const int retval(compare(target1.first, target2.first));
                // std::cout << "[0] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ") -> " << retval
                //     << std::endl;
                if (retval != 0)
                {
                    ignores_.pop_back();
                    return retval;
                }
            }

            ++i;
            ++j;
        }
        ignores_.pop_back();
        return 0;
    }

    bool operator()(const index_type& val1, const index_type& val2)
    {
        // return val1 < val2;
        ignores_.clear();
        return 0 < compare(val1, val2);
    }

    void reorder_units(
        std::vector<unsigned int>& unit_indices, const unsigned int& idx,
        unsigned int& stride)
    {
        if (unit_indices[idx] != root_.size())
        {
            return;
        }

        const UnitSpecies& usp(root_.at(idx));

        unit_indices[idx] = stride;
        ++stride;

        for (UnitSpecies::container_type::const_iterator i(usp.begin());
            i != usp.end(); ++i)
        {
            if ((*i).second.second == "" || is_wildcard((*i).second.second))
            {
                continue;
            }

            // const std::vector<species_structure::site_type>&
            //     pair((*connections_.find((*i).second.second)).second);
            const std::vector<species_structure::site_type>&
                pair(connections_[(*i).second.second]);
            const species_structure::site_type&
                tgt((pair[0].first == idx && pair[0].second == (*i).first)?
                    pair[1] : pair[0]);

            reorder_units(unit_indices, tgt.first, stride);
        }
    }

protected:

    const std::vector<UnitSpecies> root_;
    connection_container_type connections_;
    std::vector<std::pair<index_type, index_type> > ignores_;
};

// int compare_unit_species(const UnitSpecies& lhs, const UnitSpecies& rhs)
// {
//     if (lhs.name() != rhs.name())
//     {
//         return (lhs.name() < rhs.name() ? 1 : -1);
//     }
// 
//     for (UnitSpecies::container_type::const_iterator i(lhs.begin()), j(rhs.begin());
//         i != lhs.end() && j != rhs.end(); ++i, ++j)
//     {
//         if ((*i).first != (*j).first)
//         {
//             return ((*i).first < (*j).first ? 1 : -1);
//         }
//         else if ((*i).second.first != (*j).second.first)
//         {
//             return ((*i).second.first < (*j).second.first ? 1 : -1);
//         }
//         else if (((*i).second.second == "") != ((*j).second.second == ""))
//         {
//             return ((*i).second.second == "" ? 1 : -1);
//         }
//     }
// 
//     if (lhs.num_sites() != rhs.num_sites())
//     {
//         return (lhs.num_sites() < rhs.num_sites() ? 1 : -1);
//     }
//     return 0;
// }
// 
// struct less_unit_species
// {
//     bool operator()(const UnitSpecies& lhs, const UnitSpecies& rhs)
//     {
//         return compare_unit_species(lhs, rhs) > 0;
//     }
// };

Species format_species(const Species& sp)
{
    // {
    //     std::vector<Species> res;

    //     std::vector<UnitSpecies> units = sp.units();
    //     std::sort(units.begin(), units.end(), less_unit_species());

    //     species_formatter formatter(units);
    //     std::vector<species_formatter::size_type>::difference_type start = 0;
    //     while (true)
    //     {
    //         std::vector<species_formatter::size_type> permutation(units.size(), units.size());
    //         species_formatter::size_type stride = formatter.sort(permutation, start, 0);
    //         res.push_back(formatter.collect(permutation));
    //         std::vector<species_formatter::size_type>::iterator start_it = permutation.begin();
    //         std::advance(start_it, start);
    //         start = std::distance(permutation.begin(), std::find(start_it, permutation.end(), units.size()));
    //         if (start == units.size())
    //         {
    //             break;
    //         }
    //     }

    //     std::cout << "BEFORE: " << sp.serial() << std::endl;
    //     std::cout << "AFTER: ";

    //     for (std::vector<Species>::const_iterator i = res.begin(); i != res.end(); ++i)
    //     {
    //         std::cout << " " << (*i).serial();
    //     }
    //     std::cout << std::endl;
    // }

    species_structure comp(sp);
    const std::vector<UnitSpecies>::size_type num_units = comp.units().size();

    std::vector<species_structure::index_type> units;
    for (species_structure::index_type i(0); i < num_units; ++i)
    {
        units.push_back(i);
    }

    std::sort(units.begin(), units.end(), comp);

    std::vector<species_structure::index_type>
        next(num_units, num_units);
    unsigned int stride(0);
    for (species_structure::index_type i(0); i < num_units; ++i)
    {
        const species_structure::index_type idx(units[i]);
        comp.reorder_units(next, idx, stride);
    }
    for (unsigned int i(0); i < num_units; ++i)
    {
        units[next[i]] = i;
    }

    Species newsp;
    std::unordered_map<std::string, std::string> cache;
    stride = 1;
    std::stringstream ss;
    for (std::vector<species_structure::index_type>::const_iterator
        i(units.begin()); i != units.end(); ++i)
    {
        UnitSpecies usp(sp.units().at(*i));
        for (UnitSpecies::container_type::size_type j(0);
            j < static_cast<UnitSpecies::container_type::size_type>(usp.num_sites()); ++j)
        {
            UnitSpecies::container_type::value_type& site(usp.at(j));
            if (site.second.second == "" || is_wildcard(site.second.second))
            {
                continue;
            }

            std::unordered_map<std::string, std::string>::const_iterator
                it(cache.find(site.second.second));
            if (it == cache.end())
            {
                ss << stride;
                cache.insert(std::make_pair(site.second.second, ss.str()));
                site.second.second = ss.str();
                ++stride;
                ss.clear();
                ss.str("");
            }
            else
            {
                site.second.second = (*it).second;
            }
        }
        newsp.add_unit(usp);
    }
    return newsp;
}

boost::optional<rule_based_expression_matcher<UnitSpecies>::context_type>
    rule_based_expression_matcher<UnitSpecies>::match_unit_species(
        const UnitSpecies& pttrn,
        const UnitSpecies& usp,
        const rule_based_expression_matcher<UnitSpecies>::context_type& org)
{
    typedef rule_based_expression_matcher<UnitSpecies>::context_type context_type;

    context_type ctx = org;

    if (is_wildcard(pttrn.name()))
    {
        if (is_pass_wildcard(pttrn.name()))
        {
            throw NotSupported(
                "A pass wildcard '_0' is not allowed to be a name of Species.");
        }
        else if (is_named_wildcard(pttrn.name()))
        {
            context_type::variable_container_type::const_iterator
                itr(ctx.globals.find(pttrn.name()));
            if (itr == ctx.globals.end())
            {
                ctx.globals[pttrn.name()] = usp.name();
            }
            else if ((*itr).second != usp.name())
            {
                return boost::none;
            }
        }
    }
    else if (pttrn.name() != usp.name())
    {
        return boost::none;
    }

    for (UnitSpecies::container_type::const_iterator j(pttrn.begin());
        j != pttrn.end(); ++j)
    {
        if (usp.has_site((*j).first))
        {
            const UnitSpecies::site_type& site(usp.get_site((*j).first));

            if ((*j).second.first != "")
            {
                if (site.first == "")
                {
                    return boost::none;
                }
                else if (is_pass_wildcard((*j).second.first))
                {
                    throw NotSupported(
                        "A pass wildcard '_0' is not allowed to be a state.");
                }
                else if (is_unnamed_wildcard((*j).second.first))
                {
                    ; // do nothing
                }
                else if (is_named_wildcard((*j).second.first))
                {
                    context_type::variable_container_type::const_iterator
                        itr(ctx.globals.find((*j).second.first));
                    if (itr == ctx.globals.end())
                    {
                        ctx.globals[(*j).second.first] = site.first;
                    }
                    else if ((*itr).second != site.first)
                    {
                        return boost::none;
                    }
                }
                else if ((*j).second.first != site.first)
                {
                    return boost::none;
                }
            }

            if (is_pass_wildcard((*j).second.second))
            {
                ; // just skip checking
            }
            else if ((*j).second.second == "")
            {
                if (site.second != "")
                {
                    return boost::none;
                }
            }
            else
            {
                if (site.second == "")
                {
                    return boost::none;
                }
                else if (is_unnamed_wildcard((*j).second.second))
                {
                    continue;
                }
                else if (is_named_wildcard((*j).second.second))
                {
                    throw NotSupported(
                        "A named wildcard is not allowed to be a bond.");
                }

                context_type::variable_container_type::const_iterator
                    itr(ctx.locals.find((*j).second.second));
                if (itr == ctx.locals.end())
                {
                    ctx.locals[(*j).second.second] = site.second;
                }
                else if ((*itr).second != site.second)
                {
                    return boost::none;
                }

            }
        }
        else
        {
            return boost::none;
        }
    }

    return ctx;
}

/*
 * Apply a ReactionRule to the given set of reactants and return ReactionRules.
 */

struct _ReactionRuleExpressionMatcher
{
    typedef rule_based_expression_matcher<UnitSpecies>::context_type context_type;
    typedef ReactionRule::reactant_container_type reactant_container_type;

    typedef struct
    {
        std::vector<UnitSpecies> products;
        std::vector<std::vector<UnitSpecies>::size_type> correspo;
        std::vector<std::vector<UnitSpecies>::size_type> removed;
        std::vector<UnitSpecies>::size_type reserved;
    } operation_type;

    typedef struct
    {
        std::vector<UnitSpecies> units;
        std::vector<unsigned int> groups;
        unsigned int num_groups;
    } unit_group_type;
};

_ReactionRuleExpressionMatcher::operation_type compile_reaction_rule(const ReactionRule& pttrn)
{
    typedef std::vector<UnitSpecies>::size_type size_type;
    typedef std::vector<UnitSpecies>::const_iterator const_iterator;

    _ReactionRuleExpressionMatcher::operation_type res = {};
    std::vector<UnitSpecies>& products = res.products;
    std::vector<size_type>& correspo = res.correspo;
    std::vector<size_type>& removed = res.removed;

    // 1. Concatenate units of a pattern

    std::vector<UnitSpecies> reactants;
    for (ReactionRule::reactant_container_type::const_iterator
        i(pttrn.reactants().begin()); i != pttrn.reactants().end(); ++i)
    {
        std::vector<UnitSpecies> const units = (*i).units();
        reactants.reserve(reactants.size() + units.size());
        std::copy(units.begin(), units.end(), std::back_inserter(reactants));
    }

    res.reserved = reactants.size();

    int product_bond_stride = 0;
    for (ReactionRule::reactant_container_type::const_iterator
        i(pttrn.products().begin()); i != pttrn.products().end(); ++i)
    {
        product_bond_stride += concatenate_units(products, (*i), product_bond_stride);
    }

    // 2. Check correspondences between reactant and product units

    correspo.reserve(products.size());

    {
        size_type num_units(reactants.size());

        size_type idx1 = 0;
        for (const_iterator i(products.begin()); i != products.end(); ++i, ++idx1)
        {
            size_type idx2 = 0;
            for (const_iterator j(reactants.begin()); j != reactants.end(); ++j, ++idx2)
            {
                if (is_correspondent(*i, *j))
                {
                    if (correspo.size() > idx1)
                    {
                        ;  //WARN: multiple correspondence found
                        assert(false);  // never get here
                    }
                    else if (std::find(correspo.begin(), correspo.end(), idx2)
                        != correspo.end())
                    {
                        ;  //WARN: multiple correspondence skipped
                        ;  // The reactant matches is already assigned to the other product
                    }
                    else
                    {
                        correspo.push_back(idx2);
                        break;
                    }
                }
            }

            // If no reactant is found, create new one
            if (correspo.size() == idx1)
            {
                correspo.push_back(num_units);
                ++num_units;
            }
        }
    }

    // List reactants removed after reacting
    for (size_type i(0); i < reactants.size(); ++i)
    {
        if (std::find(correspo.begin(), correspo.end(), i) == correspo.end())
        {
            removed.push_back(i);
        }
    }

    return res;
}

_ReactionRuleExpressionMatcher::unit_group_type generate_units(
    const _ReactionRuleExpressionMatcher::operation_type& operations,
    const _ReactionRuleExpressionMatcher::context_type& ctx,
    const ReactionRule::reactant_container_type& target,
    const ReactionRule::policy_type& policy)
{
    typedef _ReactionRuleExpressionMatcher::context_type context_type;
    typedef _ReactionRuleExpressionMatcher::unit_group_type unit_group_type;
    typedef std::vector<UnitSpecies>::size_type size_type;

    const std::vector<UnitSpecies>& products = operations.products;
    const std::vector<size_type>& correspo = operations.correspo;
    const std::vector<size_type>& removed = operations.removed;
    const size_type& reserved = operations.reserved;

    unit_group_type res = {};
    std::vector<UnitSpecies>& units = res.units;
    std::vector<unsigned int>& groups = res.groups;
    unsigned int& num_groups = res.num_groups;

    // 3. Concatenate units given as reactants

    // const context_type ctx(context());

    int bond_stride = 0;

    //XXX: for (std::vector<reactant_container_type::size_type>::const_iterator
    //XXX:     i(permutation_.begin()); i != permutation_.end(); ++i)
    //XXX: {
    //XXX:     bond_stride += concatenate_units(units, target[*i], bond_stride);
    //XXX: }
    for (ReactionRule::reactant_container_type::const_iterator
        i(target.begin()); i != target.end(); ++i)
    {
        bond_stride += concatenate_units(units, *i, bond_stride);
    }

    // 4. Modify units

    std::unordered_map<std::string, std::string> bond_cache;

    {
        std::vector<std::pair<size_type, size_type> > priorities;
        {
            size_type next_idx = units.size();
            size_type idx = 0;
            for (std::vector<size_type>::const_iterator i(correspo.begin());
                i != correspo.end(); ++i, ++idx)
            {
                if ((*i) < reserved)
                {
                    assert(ctx.iterators.size() > (*i));
                    priorities.push_back(std::make_pair(ctx.iterators[(*i)], idx));
                }
                else
                {
                    priorities.push_back(std::make_pair(next_idx, idx));
                    ++next_idx;
                }
            }
            std::sort(priorities.begin(), priorities.end());
        }

        for (std::vector<std::pair<size_type, size_type> >::const_iterator
            itr1(priorities.begin()); itr1 != priorities.end(); ++itr1)
        {
            const UnitSpecies& op = products[(*itr1).second];
            const size_type& tgt = (*itr1).first;

            if (tgt >= units.size())
            {
                // 4-1. Create a new unit
                assert(tgt == units.size());
                // tgt = units.size();
                units.push_back(op);
                if (is_named_wildcard(op.name()))
                {
                    context_type::variable_container_type::const_iterator
                        itr(ctx.globals.find(op.name()));
                    if (itr == ctx.globals.end())
                    {
                        std::stringstream message;
                        message << "A named wildcard [" << op.name() << "] cannot be resolved.";
                        throw IllegalState(message.str());
                    }
                    else
                    {
                        units.back().set_name((*itr).second);
                    }
                }
            }
            // else
            // {
            //     tgt = ctx.iterators[idx2];
            // }

            for (UnitSpecies::container_type::const_iterator i(op.begin());
                i != op.end(); ++i)
            {
                UnitSpecies::container_type::value_type&
                    site(units[tgt].at((*i).first));

                // 4-2. Update sites
                if ((*i).second.first == "")
                {
                    ; // do nothing
                }
                else if (is_wildcard((*i).second.first))
                {
                    if ((*i).second.first.size() != 1)
                    {
                        context_type::variable_container_type::const_iterator
                            itr(ctx.globals.find((*i).second.first));
                        if (itr == ctx.globals.end())
                        {
                            std::stringstream message;
                            message << "An invalid global name [" << (*i).second.first << "] was given.";
                            throw IllegalState(message.str());
                        }
                        else
                        {
                            site.second.first = (*itr).second;
                        }
                    }
                }
                else
                {
                    site.second.first = (*i).second.first;
                }

                // 4-3. Update bonds
                if ((*i).second.second == "")
                {
                    // Just cut the bond
                    site.second.second = "";
                }
                else if (is_wildcard((*i).second.second))
                {
                    if (is_named_wildcard((*i).second.second))
                    {
                        std::stringstream message;
                        message << "A named wildcard [" << (*i).second.second << "cannot be a bond."
                            << " Only a unnamed wildcard is accepted here.";
                        throw IllegalState(message.str());
                    }

                    ; // do nothing
                }
                else
                {
                    // Make a new bond
                    std::unordered_map<std::string, std::string>::const_iterator
                        itr(bond_cache.find((*i).second.second));
                    if (itr != bond_cache.end())
                    {
                        site.second.second = (*itr).second;
                    }
                    else
                    {
                        ++bond_stride;
                        site.second.second = itos(bond_stride);
                        bond_cache.insert(std::make_pair((*i).second.second, site.second.second));
                    }
                }
            }
        }
    }

    // 5. Remove units deleted

    if (removed.size() > 0)
    {
        for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
            i(removed.begin()); i != removed.end(); ++i)
        {
            units[ctx.iterators[(*i)]] = UnitSpecies();
        }
    }

    // The following is originally in group_units()

    // 6. Check connections between bonds and units

    std::unordered_map<std::string, std::pair<std::string, unsigned int>> bondinfo;
    const unsigned int done = units.size();

    std::vector<std::vector<std::vector<UnitSpecies>::size_type> > connections;
    connections.resize(units.size());

    {
        unsigned int idx = 0;
        for (std::vector<UnitSpecies>::const_iterator i(units.begin());
            i != units.end(); ++i, ++idx)
        {
            for (UnitSpecies::container_type::const_iterator j((*i).begin());
                j != (*i).end(); ++j)
            {
                const std::string bond((*j).second.second);
                if (bond == "")
                {
                    continue;
                }
                else if (is_wildcard(bond))
                {
                    std::stringstream ss;
                    ss << "A bond in a product is a wildcard ["
                        << (*i).name() << "(" << (*j).first << "^" << bond << ")].";
                    throw IllegalState(ss.str());
                }

                std::unordered_map<std::string, std::pair<std::string, unsigned int>>::iterator
                    itr(bondinfo.find(bond));
                if (itr == bondinfo.end())
                {
                    bondinfo[bond] = std::make_pair((*j).first, idx);
                }
                else
                {
                    const unsigned int another = (*itr).second.second;
                    if (another == done)
                    {
                        std::stringstream ss;
                        ss << "A bond in a product is multiply connected ["
                            << (*i).name() << "(" << (*j).first << "^" << bond << ")].";
                        throw IllegalState(ss.str());
                    }

                    connections[idx].push_back(another);
                    connections[another].push_back(idx);
                    bondinfo[bond].second = done;  // This means the bond is already assigned.
                }
            }
        }
    }

    // const ReactionRule::policy_type& policy(pttrn.policy());

    if ((policy & ReactionRule::POLICY_STRICT)
        && !(policy & (ReactionRule::POLICY_DESTROY | ReactionRule::POLICY_IMPLICIT)))
    {
        for (std::unordered_map<std::string, std::pair<std::string, unsigned int>>::const_iterator
            i(bondinfo.begin()); i != bondinfo.end(); ++i)
        {
            if ((*i).second.second != done)
            {
                std::stringstream ss;
                ss << "A bond in a product is not resolved [" << (*i).first << "].";
                {
                    Species sp;
                    for (std::vector<UnitSpecies>::const_iterator j(units.begin());
                        j != units.end(); ++j)
                    {
                        sp.add_unit((*j));
                    }
                    ss << " " << sp.serial();
                }
                throw IllegalState(ss.str());
            }
        }
    }

    // 7. Group units based on the connections

    const unsigned int notyet = units.size();
    groups.resize(units.size(), notyet);

    for (unsigned int idx(0); idx != units.size(); ++idx)
    {
        if (units[idx] == UnitSpecies())
        {
            continue;  // A removed unit
        }

        // std::cout << idx << " connects with ";
        // for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
        //     i(connections[idx].begin()); i != connections[idx].end(); ++i)
        // {
        //     std::cout << *i;
        // }
        // std::cout << std::endl;

        num_groups = tag_units(
            groups, num_groups, idx, connections, notyet);
    }

    assert(num_groups <= notyet);

    // 8. Resolve open bonds based on the policy

    if (policy & ReactionRule::POLICY_IMPLICIT)
    {
        for (std::vector<UnitSpecies>::iterator i(units.begin());
            i != units.end(); ++i)
        {
            for (UnitSpecies::container_type::size_type j(0);
                static_cast<Integer>(j) != (*i).num_sites(); ++j)
            {
                UnitSpecies::container_type::value_type& site((*i).at(j));
                const std::string bond(site.second.second);
                if (bond != "" && !is_wildcard(bond) && bondinfo[bond].second != done)
                {
                    site.second.second = "";
                }
            }
        }
    }

    if (policy & ReactionRule::POLICY_DESTROY)
    {
        for (std::unordered_map<std::string, std::pair<std::string, unsigned int>>::const_iterator
            i(bondinfo.begin()); i != bondinfo.end(); ++i)
        {
            if ((*i).second.second != done && units[(*i).second.second] != UnitSpecies())
            {
                const unsigned int group_id = groups[(*i).second.second];

                // assert(num_groups > 0);
                // --num_groups;
                std::vector<unsigned int>::iterator j1(groups.begin());
                for (std::vector<UnitSpecies>::iterator j2(units.begin());
                    j2 != units.end(); ++j1, ++j2)
                {
                    if ((*j1) == group_id)
                    {
                        (*j2) = UnitSpecies();
                        // (*j1) = notyet;
                    }
                }
            }
        }
    }

    return res;
}

std::vector<Species> group_units(
    const std::vector<UnitSpecies>& units,
    const std::vector<unsigned int>& groups,
    const unsigned int num_groups)
{
    // 8. Divide units into Species

    std::vector<Species> products;
    products.resize(num_groups);

    // {
    //     std::vector<unsigned int>::size_type idx = 0;
    //     for (std::vector<UnitSpecies>::iterator i(units.begin());
    //         i != units.end(); ++i; ++idx)
    //     {
    //         products[groups[idx]].add_unit(*i);
    //     }
    // }

    for (unsigned int idx(0); idx != num_groups; ++idx)
    {
        std::unordered_map<std::string, std::string> new_bonds;
        unsigned int stride(1);

        for (std::vector<unsigned int>::const_iterator
            i(groups.begin()); i != groups.end(); ++i)
        {
            if (idx != *i)
            {
                continue;
            }

            UnitSpecies usp(units[std::distance(groups.begin(), i)]); //XXX: copy

            if (usp == UnitSpecies())
            {
                continue;
            }

            for (UnitSpecies::container_type::size_type j(0);
                static_cast<Integer>(j) != usp.num_sites(); ++j)
            {
                UnitSpecies::container_type::value_type&
                    site(usp.at(j));
                const std::string bond(site.second.second);
                if (bond == "" || is_wildcard(bond))
                {
                    continue;
                }

                std::unordered_map<std::string, std::string>::const_iterator
                    itr(new_bonds.find(bond));
                if (itr == new_bonds.end())
                {
                    const std::string new_bond(itos(stride));
                    ++stride;
                    new_bonds[bond] = new_bond;
                    site.second.second = new_bond;
                }
                else
                {
                    site.second.second = (*itr).second;
                }
            }

            products[idx].add_unit(usp);
        }

        // products[idx] = format_species(products[idx]);
    }

    products.erase(std::remove(products.begin(), products.end(), Species()), products.end());

    return products;
}

}; // context

std::vector<ReactionRule> generate_reaction_rules(
    const ReactionRule& pttrn,
    const ReactionRule::reactant_container_type& reactants)
{
    typedef context::_ReactionRuleExpressionMatcher::operation_type operation_type;
    typedef context::_ReactionRuleExpressionMatcher::unit_group_type unit_group_type;
    typedef context::rule_based_expression_matcher<std::vector<Species> >::context_type context_type;
    typedef std::vector<ReactionRule> return_type;

    if (pttrn.reactants().size() == 0)
    {
        return return_type(
            1, ReactionRule(reactants, pttrn.products(), pttrn.k()));  // Zeroth-order reactions
    }

    std::vector<std::vector<UnitSpecies> > candidates;
    const operation_type op = context::compile_reaction_rule(pttrn);

    return_type reactions;
    context::rule_based_expression_matcher<std::vector<Species> > matcher(pttrn.reactants());
    for (boost::optional<context_type> ctx = matcher.match(reactants); ctx; ctx = matcher.next())
    {
        const unit_group_type _res
            = context::generate_units(op, (*ctx), reactants, pttrn.policy());

        std::vector<std::vector<UnitSpecies> >::iterator
            i(std::find(candidates.begin(), candidates.end(), _res.units));
        if (i != candidates.end())
        {
            ; // (*i).set_k((*i).k() + rr.k());
        }
        else
        {
            candidates.push_back(_res.units);
            reactions.push_back(
                ReactionRule(
                    reactants,
                    context::group_units(_res.units, _res.groups, _res.num_groups),
                    pttrn.k()));
        }
    }
    return reactions;
}

} // ecell4
