#ifndef ECELL4_CONTEXT_HPP
#define ECELL4_CONTEXT_HPP

#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include <boost/array.hpp>
#include <boost/optional.hpp>


namespace ecell4
{

namespace context
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

class unit_species_comparerator
{
public:

    // typedef Species::container_type::size_type index_type;
    typedef unsigned int index_type;
    typedef std::pair<index_type, std::string> site_type;
    typedef utils::get_mapper_mf<std::string, std::vector<site_type> >::type
        connection_container_type;

public:

    unit_species_comparerator(const Species& sp)
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
                if ((*i).second.second == "" || rbex::is_wildcard((*i).second.second))
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
            if ((*i).second.second == "" || rbex::is_wildcard((*i).second.second))
            {
                continue;
            }

            // const std::vector<unit_species_comparerator::site_type>&
            //     pair((*connections_.find((*i).second.second)).second);
            const std::vector<unit_species_comparerator::site_type>&
                pair(connections_[(*i).second.second]);
            const unit_species_comparerator::site_type&
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

Species format_species(const Species& sp);

inline Species::serial_type unique_serial(const Species& sp)
{
    return format_species(sp).serial();
}

}  // context

namespace _context
{

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
            typedef utils::get_mapper_mf<std::string, std::string>::type
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
        another_ = another;
        ctx_ = ctx;

        if (const boost::optional<context_type> retval = match_unit_species(pttrn_, another_, ctx))
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
    value_type another_;
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
        : pttrn_(pttrn)
    {
        ;
    }

    boost::optional<context_type> match(const value_type& another)
    {
        context_type::variable_container_type globals;
        return match(another, globals);
    }

    boost::optional<context_type> match(
        const value_type& another, const context_type::variable_container_type& globals)
    {
        matchers_.clear();
        for (value_type::const_iterator i(pttrn_.begin());
            i != pttrn_.end(); ++i)
        {
            matchers_.push_back(submatcher_type(*i));
        }

        another_ = another;
        context_type ctx;
        ctx.globals = globals;
        return advance(0, 0, ctx);
    }

    boost::optional<context_type> next()
    {
        if (matchers_.size() == 0)
        {
            return ctx_;
        }
        return __next();
    }

    Integer count(const value_type& another)
    {
        context_type::variable_container_type globals;
        if (!match(another, globals))
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

    boost::optional<context_type> __next()
    {
        //XXX: Make sure match was already called..
        size_t src = matchers_.size();
        while (src > 0)
        {
            --src;

            submatcher_type& matcher_at_src(matchers_[src]);
            const size_t dst = ctx_.iterators[src];
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
            ctx_ = ctx;
            return ctx_;
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

    value_type pttrn_;
    value_type another_;
    std::vector<submatcher_type> matchers_;
    context_type ctx_;
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

    const context_type& context() const
    {
        return base_.context();
    }

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
        : pttrn_(pttrn)
    {
        ;
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
        return __next();
    }

    //TODO: HERE
    context_type context() const
    {
        context_type ctx;
        if (matchers_.size() == 0)
        {
            return ctx;
        }

        ctx.globals = contexts_.back().globals;

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
            const context_type& ctx_at_idx1 = contexts_[idx1];

            for (context_type::iterator_container_type::const_iterator
                j(ctx_at_idx1.iterators.begin());
                j != ctx_at_idx1.iterators.end(); ++j)
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
        permutation.reserve(pttrn_.size());
        for (size_t i(0); i != pttrn_.size(); ++i)
        {
            permutation.push_back(i);
        }
        return __match(another, permutation);
    }

    boost::optional<context_type> __match(const value_type& another, const permutation_type& permutation)
    {
        assert(0 <= pttrn_.size() < 3);
        assert(pttrn_.size() == permutation.size());

        another_ = another; //XXX: copy?
        permutation_ = permutation;

        if (pttrn_.size() != another.size())
        {
            return boost::none;
        }

        matchers_.clear();
        for (value_type::const_iterator i(pttrn_.begin()); i != pttrn_.end(); ++i)
        {
            matchers_.push_back(submatcher_type(*i));
        }

        contexts_.resize(matchers_.size());

        context_type::variable_container_type globals;
        return this->advance(0, globals);
    }

    boost::optional<context_type> __next()
    {
        //XXX: Make sure match was already called..
        size_t pos = matchers_.size();
        while (pos > 0)
        {
            --pos;

            submatcher_type& matcher_at_pos(matchers_[pos]);
            while (const boost::optional<context_type> res = matcher_at_pos.next())
            {
                contexts_[pos] = (*res);
                if (const boost::optional<context_type> newres = this->advance(pos + 1, (*res).globals))
                {
                    return newres;
                }
            }
        }

        return boost::none;
    }

    boost::optional<context_type> advance(const size_t pos, const context_type::variable_container_type& globals)
    {
        if (pos == matchers_.size())
        {
            return this->context();
        }

        submatcher_type& matcher_at_pos(matchers_[pos]);
        const value_type::value_type& another_at_pos(another_[pos]);
        for (boost::optional<context_type> res = matcher_at_pos.match(another_at_pos, globals);
             res; res = matcher_at_pos.next())
        {
            contexts_[pos] = (*res);
            if (const boost::optional<context_type> newres = this->advance(pos + 1, (*res).globals))
            {
                return newres;
            }
        }
        return boost::none;
    }

protected:

    value_type pttrn_;
    value_type another_;
    permutation_type permutation_;

    std::vector<submatcher_type> matchers_;
    std::vector<context_type> contexts_;
};

} // _context

/**
 * New interfaces for the rule-based modeling
 */

inline Species format_species(const Species& sp)
{
    return context::format_species(sp);
}

struct SpeciesExpressionMatcher
{
    _context::rule_based_expression_matcher<Species> pttrn;

    SpeciesExpressionMatcher(const Species& pttrn)
        : pttrn(pttrn)
    {
        ;
    }

    bool match(const Species& sp)
    {
        return pttrn.match(sp);
    }

    bool next()
    {
        return pttrn.next();
    }

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
