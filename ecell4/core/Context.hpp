#ifndef ECELL4_CONTEXT_HPP
#define ECELL4_CONTEXT_HPP

#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include <boost/array.hpp>


namespace ecell4
{

namespace rbex
{

inline bool is_wildcard(const std::string& name)
{
    return (name.size() > 0 && name[0] == '_');
}

inline bool is_unnamed_wildcard(const std::string& name)
{
    return name == "_";
}

inline bool is_pass_wildcard(const std::string& name)
{
    return name == "_0";
}

inline bool is_named_wildcard(const std::string& name)
{
    return (name.size() > 1 && name[0] == '_' && !is_pass_wildcard(name));
}

} // rbex

class MatchObject
{
public:

    typedef struct
    {
        typedef std::vector<Species::container_type::difference_type>
            iterator_container_type;
        typedef utils::get_mapper_mf<std::string, std::string>::type
            variable_container_type;

        iterator_container_type iterators;
        variable_container_type locals;
        variable_container_type globals;
    } context_type;

public:

    MatchObject(const UnitSpecies& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    virtual ~MatchObject()
    {
        ;
    }

    std::pair<bool, context_type> match(
        const Species& sp, const context_type& ctx)
    {
        // target_ = sp;
        target_ = sp.units();
        itr_ = target_.begin();
        ctx_ = ctx;
        return next();
    }

    std::pair<bool, context_type> match(
        const std::vector<UnitSpecies>& target, const context_type& ctx)
    {
        target_ = target;
        itr_ = target_.begin();
        ctx_ = ctx;
        return next();
    }

    std::pair<bool, context_type> next();

protected:

    UnitSpecies pttrn_;
    // Species target_;
    std::vector<UnitSpecies> target_;
    Species::container_type::const_iterator itr_;
    context_type ctx_;
};

std::pair<bool, MatchObject::context_type>
uspmatch(const UnitSpecies& pttrn, const UnitSpecies& sp,
    const MatchObject::context_type& org);
bool __spmatch(
    Species::container_type::const_iterator itr,
    const Species::container_type::const_iterator& end,
    const Species& sp, const MatchObject::context_type& ctx);
bool spmatch(const Species& pttrn, const Species& sp);
Integer count_spmatches(const Species& pttrn, const Species& sp);
Integer count_spmatches(
    const Species& pttrn, const Species& sp,
    const MatchObject::context_type::variable_container_type& globals);

// ReactionRule create_reaction_rule_formatted(
//     const ReactionRule::reactant_container_type& reactants,
//     const ReactionRule::product_container_type& products, const Real k);

// inline ReactionRule create_reaction_rule_formatted(const ReactionRule& rr)
// {
//     return create_reaction_rule_formatted(rr.reactants(), rr.products(), rr.k());
// }

class SpeciesExpressionMatcher
{
public:

    typedef MatchObject::context_type context_type;

public:

    SpeciesExpressionMatcher(const Species& pttrn)
        : pttrn_(pttrn.units())
    {
        ;
    }

    virtual ~SpeciesExpressionMatcher()
    {
        ;
    }

    bool match(const Species& sp)
    {
        context_type::variable_container_type globals;
        return match(sp, globals);
    }

    bool match(
        const Species& sp, const context_type::variable_container_type& globals)
    {
        matches_.clear();
        for (Species::container_type::const_iterator i(pttrn_.begin());
            i != pttrn_.end(); ++i)
        {
            matches_.push_back(MatchObject(*i));
        }

        // target_ = sp;
        target_ = sp.units();
        itr_ = matches_.begin();
        context_type ctx;
        ctx.globals = globals;
        return __match(ctx);
    }

    bool __match(const context_type& ctx)
    {
        if (itr_ == matches_.end())
        {
            ctx_ = ctx;
            return true;
        }

        std::pair<bool, context_type> retval((*itr_).match(target_, ctx));
        while (retval.first)
        {
            ++itr_;
            const bool succeeded(__match(retval.second));
            if (succeeded)
            {
                return true;
            }
            --itr_;
            retval = (*itr_).next();
        }
        return false;
    }

    bool next()
    {
        if (itr_ != matches_.end())
        {
            return false;
        }
        else if (matches_.size() == 0)
        {
            return true;
        }

        do
        {
            --itr_;
            std::pair<bool, context_type> retval((*itr_).next());
            while (retval.first)
            {
                ++itr_;
                const bool succeeded(__match(retval.second));
                if (succeeded)
                {
                    return true;
                }
                --itr_;
                retval = (*itr_).next();
            }
        }
        while (itr_ != matches_.begin());
        return false;
    }

    Integer count(const Species& sp)
    {
        context_type::variable_container_type globals;
        if (!match(sp, globals))
        {
            return 0;
        }
        Integer n(1);
        while (next())
        {
            ++n;
        }
        return n;
    }

    const context_type& context() const
    {
        return ctx_;
    }

protected:

    // Species pttrn_;
    // Species target_;
    std::vector<UnitSpecies> pttrn_;
    std::vector<UnitSpecies> target_;
    std::vector<MatchObject> matches_;
    std::vector<MatchObject>::iterator itr_;
    context_type ctx_;
};

class ReactionRuleExpressionMatcher
{
public:

    typedef MatchObject::context_type context_type;
    typedef ReactionRule::reactant_container_type reactant_container_type;

public:

    ReactionRuleExpressionMatcher(const ReactionRule& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    virtual ~ReactionRuleExpressionMatcher()
    {
        ;
    }

    bool match(const Species& sp)
    {
        reactant_container_type reactants;
        reactants.push_back(sp);
        permutation_.clear();
        permutation_.push_back(0);
        return __match(reactants);
    }

    bool match(const Species& sp1, const Species& sp2)
    {
        reactant_container_type reactants;
        reactants.push_back(sp1);
        reactants.push_back(sp2);
        permutation_.clear();
        permutation_.push_back(0);
        permutation_.push_back(1);
        return __match(reactants);
    }

    bool match_reversed(const Species& sp1, const Species& sp2)
    {
        reactant_container_type reactants;
        reactants.push_back(sp2);
        reactants.push_back(sp1);
        permutation_.clear();
        permutation_.push_back(1);
        permutation_.push_back(0);
        return __match(reactants);
    }

    bool match(const reactant_container_type& reactants,
               const std::vector<reactant_container_type::size_type>& permutation)
    {
        permutation_ = permutation;
        return __match(reactants);
    }

    bool match(const reactant_container_type& reactants)
    {
        permutation_.clear();
        permutation_.reserve(reactants.size());
        for (std::vector<reactant_container_type::size_type>::size_type i(0);
            i != reactants.size(); ++i)
        {
            permutation_.push_back(i);
        }
        return __match(reactants);
    }

    bool __match(const reactant_container_type& reactants)
    {
        if (pttrn_.reactants().size() != reactants.size())
        {
            return false;
        }

        matchers_.clear();
        for (reactant_container_type::const_iterator
            i(pttrn_.reactants().begin()); i != pttrn_.reactants().end(); ++i)
        {
            matchers_.push_back(SpeciesExpressionMatcher(*i));
        }

        target_ = reactants; //XXX: copy?
        itr_ = matchers_.begin();
        context_type::variable_container_type globals;
        return __submatch(globals);
    }

    bool __submatch(const context_type::variable_container_type& globals)
    {
        if (itr_ == matchers_.end())
        {
            return true;
        }

        bool retval((*itr_).match(
            target_[std::distance(matchers_.begin(), itr_)],
            globals));
        while (retval)
        {
            const context_type::variable_container_type&
                globals_prev((*itr_).context().globals);
            ++itr_;
            const bool succeeded(__submatch(globals_prev));
            if (succeeded)
            {
                return true;
            }
            --itr_;
            retval = (*itr_).next();
        }
        return false;
    }

    bool next()
    {
        if (itr_ != matchers_.end() || pttrn_.reactants().size() == 0)
        {
            return false;
        }
        else if (matchers_.size() == 0)
        {
            return true;
        }

        do
        {
            --itr_;
            bool retval((*itr_).next());
            while (retval)
            {
                const context_type::variable_container_type&
                    globals_prev((*itr_).context().globals);
                ++itr_;
                const bool succeeded(__submatch(globals_prev));
                if (succeeded)
                {
                    return true;
                }
                --itr_;
                retval = (*itr_).next();
            }
        }
        while (itr_ != matchers_.begin());
        return false;
    }

    std::pair<bool, context_type> __match(
        const context_type::variable_container_type& globals,
        reactant_container_type::const_iterator i,
        reactant_container_type::const_iterator j)
    {
        SpeciesExpressionMatcher m(*i);
        if (!m.match(*j, globals))
        {
            return std::make_pair(false, context_type());
        }

        ++i;
        ++j;
        if (i == pttrn_.reactants().end() || j == target_.end())
        {
            return std::make_pair(true, m.context());
        }

        do
        {
            if (__match(m.context().globals, i, j).first)
            {
                return std::make_pair(true, m.context());
            }
        } while (m.next());
        return std::make_pair(false, context_type());
    }

    context_type context() const
    {
        context_type ctx;
        if (matchers_.size() == 0)
        {
            return ctx;
        }

        ctx.globals = matchers_.back().context().globals;

        std::vector<unsigned int> strides;
        strides.reserve(target_.size());
        {
            unsigned int stride = 0;
            for (std::vector<Species>::const_iterator
                i(target_.begin()); i != target_.end(); ++i)
            {
                strides.push_back(stride);
                stride += (*i).units().size();
            }
        }

        for (std::vector<SpeciesExpressionMatcher>::const_iterator
            i(matchers_.begin()); i != matchers_.end(); ++i)
        {
            const unsigned int idx1 = std::distance(matchers_.begin(), i);  // a position in matcher_
            const unsigned int idx2 = permutation_[idx1];  // a position in reactants
            const unsigned int stride = strides[idx2];

            for (context_type::iterator_container_type::const_iterator
                j((*i).context().iterators.begin());
                j != (*i).context().iterators.end(); ++j)
            {
                // const unsigned int idx3 = std::distance((*i).context().iterators.begin(), j);  // a position in context.iterators
                const unsigned int idx4 = (*j);  // a position in units of a Species

                ctx.iterators.push_back(idx4 + stride);
            }
        }

        // Species::container_type::difference_type stride(0);
        // for (std::vector<SpeciesExpressionMatcher>::const_iterator
        //     i(matchers_.begin()); i != matchers_.end(); ++i)
        // {
        //     for (context_type::iterator_container_type::const_iterator
        //         j((*i).context().iterators.begin());
        //         j != (*i).context().iterators.end(); ++j)
        //     {
        //         ctx.iterators.push_back((*j) + stride);
        //     }
        //     stride += target_[std::distance(matchers_.begin(), i)].units().size();
        // }
        //XXX: Species::container_type::difference_type totstride(0);
        //XXX: std::vector<Species::container_type::difference_type> strides(matchers_.size());
        //XXX: for (std::vector<reactant_container_type::size_type>::const_iterator
        //XXX:     i(permutation_.begin()); i != permutation_.end(); ++i)
        //XXX: {
        //XXX:     strides[(*i)] = totstride;
        //XXX:     totstride += target_[(*i)].units().size();
        //XXX: }

        //XXX: for (std::vector<SpeciesExpressionMatcher>::const_iterator
        //XXX:     i(matchers_.begin()); i != matchers_.end(); ++i)
        //XXX: {
        //XXX:     const Species::container_type::difference_type stride
        //XXX:         = strides[std::distance(matchers_.begin(), i)];
        //XXX:     for (context_type::iterator_container_type::const_iterator
        //XXX:         j((*i).context().iterators.begin());
        //XXX:         j != (*i).context().iterators.end(); ++j)
        //XXX:     {
        //XXX:         ctx.iterators.push_back((*j) + stride);
        //XXX:     }
        //XXX: }

        return ctx;
    }

    std::vector<Species> generate();

    typedef struct
    {
        std::vector<UnitSpecies> products;
        std::vector<std::vector<UnitSpecies>::size_type> correspo;
        std::vector<std::vector<UnitSpecies>::size_type> removed;
        std::vector<UnitSpecies>::size_type reserved;
    } operation_type;

    operation_type compile();

    typedef struct
    {
        std::vector<UnitSpecies> units;
        std::vector<unsigned int> groups;
        unsigned int num_groups;
    } unit_group_type;

    unit_group_type genunits(const operation_type& op);
    std::vector<ReactionRule> gen(const ReactionRule::reactant_container_type& reactants);

    const reactant_container_type& reactants() const
    {
        return target_;
    }

protected:

    const ReactionRule pttrn_;
    reactant_container_type target_;
    std::vector<reactant_container_type::size_type> permutation_;
    std::vector<SpeciesExpressionMatcher> matchers_;
    std::vector<SpeciesExpressionMatcher>::iterator itr_;
};

std::vector<Species> group_units(
    const std::vector<UnitSpecies>& units,
    const std::vector<unsigned int>& groups, const unsigned int num_groups);

} // ecell4

#endif /* ECELL4_CONTEXT_HPP */
