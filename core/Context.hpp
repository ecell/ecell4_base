#ifndef __ECELL4_CONTEXT_HPP
#define __ECELL4_CONTEXT_HPP

#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include <boost/array.hpp>


namespace ecell4
{

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
        target_ = sp;
        itr_ = target_.begin();
        ctx_ = ctx;
        return next();
    }

    std::pair<bool, context_type> next();

protected:

    UnitSpecies pttrn_;
    Species target_;
    Species::container_type::const_iterator itr_;
    context_type ctx_;
};

bool is_wildcard(const std::string& name);
bool is_unnamed_wildcard(const std::string& name);
bool is_named_wildcard(const std::string& name);

std::pair<bool, MatchObject::context_type>
uspmatch(const UnitSpecies& pttrn, const UnitSpecies& sp,
    const MatchObject::context_type& org);
bool __spmatch(
    Species::container_type::const_iterator itr,
    const Species::container_type::const_iterator& end,
    const Species& sp, const MatchObject::context_type& ctx);
bool spmatch(const Species& pttrn, const Species& sp);
bool rrmatch(boost::array<Species, 2> pttrn, boost::array<Species, 2> reactants);
Integer count_spmatches(const Species& pttrn, const Species& sp);
Integer count_spmatches(
    const Species& pttrn, const Species& sp,
    const MatchObject::context_type::variable_container_type& globals);

class SpeciesExpressionMatcher
{
public:

    typedef MatchObject::context_type context_type;

public:

    SpeciesExpressionMatcher(const Species& pttrn)
        : pttrn_(pttrn)
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

        target_ = sp;
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

    const context_type& context() const
    {
        return ctx_;
    }

protected:

    const Species pttrn_;
    Species target_;
    std::vector<MatchObject> matches_;
    std::vector<MatchObject>::iterator itr_;
    context_type ctx_;
};

} // ecell4

#endif /* __ECELL4_CONTEXT_HPP */
