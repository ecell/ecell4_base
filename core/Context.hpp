#ifndef __ECELL4_CONTEXT_HPP
#define __ECELL4_CONTEXT_HPP

#include "get_mapper_mf.hpp"
#include "Species.hpp"


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

    const UnitSpecies pttrn_;
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
bool spmatch(const Species& pttrn, const Species& sp);

} // ecell4

#endif /* __ECELL4_CONTEXT_HPP */
