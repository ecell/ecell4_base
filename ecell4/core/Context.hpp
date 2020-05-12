#ifndef ECELL4_CONTEXT_HPP
#define ECELL4_CONTEXT_HPP

#include "UnitSpecies.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include <boost/optional.hpp>
#include <unordered_map>

namespace ecell4
{

namespace context
{

inline bool is_empty(const std::string& name)
{
    return name == "";
}

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

Species format_species(const Species& sp);

inline Species::serial_type unique_serial(const Species& sp)
{
    return format_species(sp).serial();
}

template <typename T>
class rule_based_expression_matcher {};

template <>
class rule_based_expression_matcher<UnitSpecies>
{
public:

    typedef UnitSpecies value_type;

    typedef struct
        {
            typedef std::vector<Species::container_type::difference_type>
                iterator_container_type;
            typedef std::unordered_map<std::string, std::string>
                variable_container_type;

            iterator_container_type iterators;
            variable_container_type locals;
            variable_container_type globals;
        } context_type;

public:

    rule_based_expression_matcher(const value_type& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    boost::optional<context_type> match(const value_type& another)
    {
        return match(another, context_type());
    }

    boost::optional<context_type> match(const value_type& another, const context_type& ctx)
    {
        // another_ = another;
        ctx_ = ctx;

        if (const boost::optional<context_type> retval = match_unit_species(pttrn_, another, ctx))
        {
            return retval;
        }

        return boost::none;
    }

    boost::optional<context_type> next()
    {
        return boost::none;
    }

    const context_type& context() const
    {
        return ctx_;
    }

protected:

    static boost::optional<context_type> match_unit_species(
        const UnitSpecies& pttrn, const UnitSpecies& usp, const context_type& org);

protected:

    value_type pttrn_;
    // value_type another_;
    context_type ctx_;
};

template <>
class rule_based_expression_matcher<std::vector<UnitSpecies> >
{
public:

    typedef std::vector<UnitSpecies> value_type;
    typedef rule_based_expression_matcher<UnitSpecies> submatcher_type;

    typedef submatcher_type::context_type context_type;

public:

    rule_based_expression_matcher(const value_type& pttrn)
    {
        matchers_.clear();
        for (value_type::const_iterator i(pttrn.begin());
            i != pttrn.end(); ++i)
        {
            matchers_.push_back(submatcher_type(*i));
        }
    }

    boost::optional<context_type> match(const value_type& another)
    {
        context_type::variable_container_type globals;
        return match(another, globals);
    }

    boost::optional<context_type> match(
        const value_type& another, const context_type::variable_container_type& globals)
    {
        another_ = another;
        context_type ctx;
        ctx.globals = globals;
        return advance(0, 0, ctx);
    }

    boost::optional<context_type> next()
    {
        if (matchers_.size() == 0)
        {
            return context_type();
        }
        return __next();
    }

    Integer count(const value_type& another)
    {
        context_type::variable_container_type globals;

        std::vector<context_type::iterator_container_type> results;
        boost::optional<context_type> ctx = match(another, globals);

        if (!ctx)
        {
            return 0;
        }
        results.push_back(ctx.get().iterators);
        std::sort(results.back().begin(), results.back().end());

        while((ctx = next()))
        {
            results.push_back(ctx.get().iterators);
            std::sort(results.back().begin(), results.back().end());
        }

        return static_cast<Integer>(std::distance(results.begin(), std::unique(results.begin(), results.end())));

        // if (!match(another, globals))
        // {
        //     return 0;
        // }

        // Integer n(1);
        // while (next())
        // {
        //     ++n;
        // }
        // return n;
    }

    // const context_type& context() const
    // {
    //     return ctx_;
    // }

protected:

    boost::optional<context_type> __next()
    {
        //XXX: Make sure match was already called..
        size_t src = matchers_.size();
        while (src > 0)
        {
            --src;

            submatcher_type& matcher_at_src(matchers_[src]);
            const size_t dst = iterators_[src];
            if (const boost::optional<context_type> res = this->advance(src, dst + 1, matcher_at_src.context()))
            {
                return res;
            }
        }
        return boost::none;
    }

    boost::optional<context_type> advance(const size_t src, const size_t dst, const context_type& ctx)
    {
        if (src == matchers_.size())
        {
            iterators_ = ctx.iterators;
            return ctx;
        }
        else if (dst == another_.size())
        {
            return boost::none;
        }
        else if (std::find(ctx.iterators.begin(), ctx.iterators.end(), dst)
                 != ctx.iterators.end())
        {
            return this->advance(src, dst + 1, ctx);
        }

        submatcher_type& matcher_at_src(matchers_[src]);
        boost::optional<context_type> res = matcher_at_src.match(another_[dst], ctx);
        if (res)
        {
            (*res).iterators.push_back(dst);

            if (const boost::optional<context_type> newres = this->advance(src + 1, 0, (*res)))
            {
                return newres;
            }
        }
        return this->advance(src, dst + 1, ctx);
    }

protected:

    value_type another_;
    std::vector<submatcher_type> matchers_;
    context_type::iterator_container_type iterators_;
};

template <>
class rule_based_expression_matcher<Species>
{
public:

    typedef Species value_type;
    typedef rule_based_expression_matcher<std::vector<UnitSpecies> > base_type;
    typedef base_type::context_type context_type;

    rule_based_expression_matcher(const value_type& pttrn)
        : base_(pttrn.units())
    {
        ;
    }

    boost::optional<context_type> match(const value_type& another)
    {
        return base_.match(another.units());
    }

    boost::optional<context_type> match(
        const value_type& another, const context_type::variable_container_type& globals)
    {
        return base_.match(another.units(), globals);
    }

    boost::optional<context_type> next()
    {
        return base_.next();
    }

    size_t count(const value_type& another)
    {
        return base_.count(another.units());
    }

    // const context_type& context() const
    // {
    //     return base_.context();
    // }

protected:

    base_type base_;
};

template <>
class rule_based_expression_matcher<std::vector<Species> >
{
public:

    typedef std::vector<Species> value_type;
    typedef rule_based_expression_matcher<Species> submatcher_type;
    typedef submatcher_type::context_type context_type;

public:

    rule_based_expression_matcher(const value_type& pttrn)
    {
        matchers_.clear();
        for (value_type::const_iterator i(pttrn.begin()); i != pttrn.end(); ++i)
        {
            matchers_.push_back(submatcher_type(*i));
        }

        iterators_.resize(matchers_.size());
        // contexts_.resize(matchers_.size());
    }

    boost::optional<context_type> match(const value_type& another)
    {
        return __match(another);
    }

    boost::optional<context_type> match(
        const value_type& another, const context_type::variable_container_type& globals)
    {
        //XXX: Not implemented yet
        return __match(another);
    }

    boost::optional<context_type> next()
    {
        //XXX: Make sure match was already called..
        size_t pos = matchers_.size();
        while (pos > 0)
        {
            --pos;

            submatcher_type& matcher_at_pos(matchers_[pos]);
            while (const boost::optional<context_type> res = matcher_at_pos.next())
            {
                // contexts_[pos] = (*res);
                iterators_[pos] = (*res).iterators;
                if (const boost::optional<context_type> newres = this->advance(pos + 1, (*res).globals))
                {
                    return newres;
                }
            }
        }
        return boost::none;
    }

    //TODO: HERE
    context_type context() const
    {
        context_type ctx;
        if (matchers_.size() == 0)
        {
            return ctx;
        }

        // ctx.globals = contexts_.back().globals;
        ctx.globals = globals_;

        std::vector<unsigned int> strides;
        strides.reserve(another_.size());
        {
            unsigned int stride = 0;
            for (value_type::const_iterator i(another_.begin()); i != another_.end(); ++i)
            {
                strides.push_back(stride);
                stride += (*i).units().size();
            }
        }

        for (std::vector<submatcher_type>::const_iterator
            i(matchers_.begin()); i != matchers_.end(); ++i)
        {
            const unsigned int idx1 = std::distance(matchers_.begin(), i);  // a position in matcher_
            const unsigned int idx2 = permutation_[idx1];  // a position in reactants
            const unsigned int stride = strides[idx2];
            // const context_type& ctx_at_idx1 = contexts_[idx1];
            const context_type::iterator_container_type& it = iterators_[idx1];

            // for (context_type::iterator_container_type::const_iterator
            //     j(ctx_at_idx1.iterators.begin());
            //     j != ctx_at_idx1.iterators.end(); ++j)
            for (context_type::iterator_container_type::const_iterator
                j(it.begin()); j != it.end(); ++j)
            {
                // // a position in context.iterators
                // const unsigned int idx3 = std::distance((*i).context().iterators.begin(), j);
                const unsigned int idx4 = (*j);  // a position in units of a Species

                ctx.iterators.push_back(idx4 + stride);
            }
        }

        return ctx;
    }
    //TODO: THERE

protected:

    typedef std::vector<value_type::size_type> permutation_type;

    boost::optional<context_type> __match(const value_type& another)
    {
        permutation_type permutation;
        permutation.reserve(matchers_.size());
        for (size_t i(0); i != matchers_.size(); ++i)
        {
            permutation.push_back(i);
        }
        return __match(another, permutation);
    }

    boost::optional<context_type> __match(const value_type& another, const permutation_type& permutation)
    {
        assert(matchers_.size() < 3);
        assert(matchers_.size() == permutation.size());

        another_ = another; //XXX: copy?
        permutation_ = permutation;

        if (matchers_.size() != another.size())
        {
            return boost::none;
        }

        context_type::variable_container_type globals;
        return this->advance(0, globals);
    }

    boost::optional<context_type> advance(const size_t pos, const context_type::variable_container_type& globals)
    {
        if (pos == matchers_.size())
        {
            globals_ = globals;
            return this->context();
        }

        submatcher_type& matcher_at_pos(matchers_[pos]);
        const value_type::value_type& another_at_pos(another_[pos]);
        for (boost::optional<context_type> res = matcher_at_pos.match(another_at_pos, globals);
             res; res = matcher_at_pos.next())
        {
            // contexts_[pos] = (*res);
            iterators_[pos] = (*res).iterators;
            if (const boost::optional<context_type> newres = this->advance(pos + 1, (*res).globals))
            {
                return newres;
            }
        }
        return boost::none;
    }

protected:

    // value_type pttrn_;
    value_type another_;
    permutation_type permutation_;

    std::vector<submatcher_type> matchers_;

    // std::vector<context_type> contexts_;
    std::vector<context_type::iterator_container_type> iterators_;
    context_type::variable_container_type globals_;
};

} // context

/**
 * New interfaces for the rule-based modeling
 */

inline Integer count_species_matches(const Species& pttrn, const Species& sp)
{
    return static_cast<Integer>(
        context::rule_based_expression_matcher<Species>(pttrn).count(sp));
}

inline Species format_species(const Species& sp)
{
    return context::format_species(sp);
}

struct SpeciesExpressionMatcher
{
    context::rule_based_expression_matcher<Species> pttrn;

    SpeciesExpressionMatcher(const Species& pttrn)
        : pttrn(pttrn)
    {
        ;
    }

    bool match(const Species& sp)
    {
        return (pttrn.match(sp) ? true : false);
    }

    // bool next()
    // {
    //     return (pttrn.next() ? true : false);
    // }

    size_t count(const Species& sp)
    {
        return pttrn.count(sp);
    }
};

std::vector<ReactionRule> generate_reaction_rules(
    const ReactionRule& pttrn,
    const ReactionRule::reactant_container_type& reactants);

} // ecell4

#endif /* ECELL4_CONTEXT_HPP */
