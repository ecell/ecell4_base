#include "Context.hpp"


namespace ecell4
{

std::pair<bool, MatchObject::context_type> uspmatch(
    const UnitSpecies& pttrn, const UnitSpecies& usp,
    MatchObject::context_type& ctx)
{
    if (pttrn.name() != "_" && pttrn.name() != usp.name())
    {
        return std::make_pair(false, ctx);
    }

    for (UnitSpecies::container_type::const_iterator j(pttrn.begin());
        j != pttrn.end(); ++j)
    {
        if (usp.has_site((*j).first))
        {
            const UnitSpecies::site_type& site(usp.get_site((*j).first));
            if ((*j).second.first != "" && (*j).second.first != site.first)
            {
                return std::make_pair(false, ctx);
            }
            else if ((*j).second.second == "" && site.second != "")
            {
                return std::make_pair(false, ctx);
            }
            else if ((*j).second.second != "")
            {
                if (site.second == "")
                {
                    return std::make_pair(false, ctx);
                }
                else if ((*j).second.second == "_")
                {
                    continue;
                }

                MatchObject::context_type::variable_container_type::const_iterator
                    itr(ctx.locals.find((*j).second.second));
                if (itr == ctx.locals.end())
                {
                    ctx.locals.insert(
                        std::make_pair((*j).second.second, site.second));
                    // ctx.locals[(*j).second.second] = site.second;
                }
                else if ((*itr).second != site.second)
                {
                    return std::make_pair(false, ctx);
                }
            }
        }
        else
        {
            return std::make_pair(false, ctx);
        }
    }
    return std::make_pair(true, ctx);
}

bool __spmatch(
    Species::container_type::const_iterator itr,
    const Species::container_type::const_iterator& end,
    const Species& sp, MatchObject::context_type ctx)
{
    if (itr == end)
    {
        for (MatchObject::context_type::iterator_container_type::const_iterator
            i(ctx.iterators.begin()); i != ctx.iterators.end(); ++i)
            std::cout << *i << " ";
        std::cout << std::endl;
        return true;
    }

    MatchObject obj(*itr);
    ++itr;

    std::pair<bool, MatchObject::context_type> retval(obj.match(sp, ctx));
    while (retval.first)
    {
        if (__spmatch(itr, end, sp, retval.second))
        {
            return true;
        }
        retval = obj.next();
    }
    return false;
}

bool spmatch(const Species& pttrn, const Species& sp)
{
    MatchObject::context_type ctx;
    return __spmatch(pttrn.begin(), pttrn.end(), sp, ctx);
}

std::pair<bool, MatchObject::context_type> MatchObject::next()
{
    for (; itr_ != target_.end(); ++itr_)
    {
        const Species::container_type::difference_type
            pos(distance(target_.begin(), itr_));
        if (std::find(ctx_.iterators.begin(), ctx_.iterators.end(), pos)
            != ctx_.iterators.end())
        {
            continue;
        }

        const UnitSpecies& usp(*itr_);
        MatchObject::context_type ctx(ctx_);
        ctx.iterators.push_back(pos);
        const std::pair<bool, MatchObject::context_type>
            retval(uspmatch(pttrn_, usp, ctx));
        if (retval.first)
        {
            ++itr_;
            return retval;
        }
    }
    return std::make_pair(false, MatchObject::context_type());
}

} // ecell4
