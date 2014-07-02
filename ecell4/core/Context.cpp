#include "Context.hpp"


namespace ecell4
{

bool is_wildcard(const std::string& name)
{
    return (name.size() > 0 && name[0] == '_');
}

bool is_unnamed_wildcard(const std::string& name)
{
    return name == "_";
}

bool is_named_wildcard(const std::string& name)
{
    return (name.size() > 1 && name[0] == '_');
}

bool is_freebond(const std::string& name)
{
    return name == "_free";
}

std::pair<bool, MatchObject::context_type> uspmatch(
    const UnitSpecies& pttrn, const UnitSpecies& usp,
    const MatchObject::context_type& org)
{
    std::pair<bool, MatchObject::context_type>
        retval(std::make_pair(false, org));
    MatchObject::context_type& ctx(retval.second);

    if (is_wildcard(pttrn.name()))
    {
        if (is_named_wildcard(pttrn.name()))
        {
            MatchObject::context_type::variable_container_type::const_iterator
                itr(ctx.globals.find(pttrn.name()));
            if (itr == ctx.globals.end())
            {
                ctx.globals[pttrn.name()] = usp.name();
            }
            else if ((*itr).second != usp.name())
            {
                return retval;
            }
        }
    }
    else if (pttrn.name() != usp.name())
    {
        return retval;
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
                    return retval;
                }
                else if (is_unnamed_wildcard((*j).second.first))
                {
                    ; // do nothing
                }
                else if (is_named_wildcard((*j).second.first))
                {
                    MatchObject::context_type::variable_container_type::const_iterator
                        itr(ctx.globals.find((*j).second.first));
                    if (itr == ctx.globals.end())
                    {
                        ctx.globals[(*j).second.first] = site.first;
                    }
                    else if ((*itr).second != site.first)
                    {
                        return retval;
                    }
                }
                else if ((*j).second.first != site.first)
                {
                    return retval;
                }
            }

            if (is_freebond((*j).second.second))
            {
                ; // just skip checking
            }
            else if ((*j).second.second == "")
            {
                if (site.second != "")
                {
                    return retval;
                }
            }
            else
            {
                if (site.second == "")
                {
                    return retval;
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

                MatchObject::context_type::variable_container_type::const_iterator
                    itr(ctx.locals.find((*j).second.second));
                if (itr == ctx.locals.end())
                {
                    ctx.locals[(*j).second.second] = site.second;
                }
                else if ((*itr).second != site.second)
                {
                    return retval;
                }

            }
        }
        else
        {
            return retval;
        }
    }

    retval.first = true;
    return retval;
}

bool __spmatch(
    Species::container_type::const_iterator itr,
    const Species::container_type::const_iterator& end,
    const Species& sp, const MatchObject::context_type& ctx)
{
    if (itr == end)
    {
        // for (MatchObject::context_type::iterator_container_type::const_iterator
        //     i(ctx.iterators.begin()); i != ctx.iterators.end(); ++i)
        //     std::cout << *i << " ";
        // std::cout << std::endl;
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
    SpeciesExpressionMatcher sexp(pttrn);
    return sexp.match(sp);
    // MatchObject::context_type ctx;
    // return __spmatch(pttrn.begin(), pttrn.end(), sp, ctx);
}

Integer count_spmatches(const Species& pttrn, const Species& sp)
{
    MatchObject::context_type::variable_container_type globals;
    return count_spmatches(pttrn, sp, globals);
}

Integer count_spmatches(const Species& pttrn, const Species& sp,
    const MatchObject::context_type::variable_container_type& globals)
{
    SpeciesExpressionMatcher sexp(pttrn);
    if (!sexp.match(sp, globals))
    {
        return 0;
    }
    Integer n(1);
    while (sexp.next())
    {
        ++n;
    }
    return n;
}

std::pair<bool, MatchObject::context_type> __rrmatch(
    const ReactionRule& rr,
    const ReactionRule::reactant_container_type& reactants,
    const MatchObject::context_type::variable_container_type& globals,
    ReactionRule::reactant_container_type::const_iterator i,
    ReactionRule::reactant_container_type::const_iterator j)
{
    SpeciesExpressionMatcher m(*i);
    if (!m.match(*j, globals))
    {
        return std::make_pair(false, MatchObject::context_type());
    }

    ++i;
    ++j;
    if (i == rr.reactants().end() || j == reactants.end())
    {
        return std::make_pair(true, m.context());
    }

    do
    {
        if (__rrmatch(rr, reactants, m.context().globals, i, j).first)
        {
            return std::make_pair(true, m.context());
        }
    } while (m.next());
    return std::make_pair(false, MatchObject::context_type());
}

bool rrmatch(const ReactionRule& rr,
    const ReactionRule::reactant_container_type& reactants)
{
    if (rr.reactants().size() != reactants.size())
    {
        return false;
    }

    ReactionRule::reactant_container_type::const_iterator
        i(rr.reactants().begin()), j(reactants.begin());
    MatchObject::context_type::variable_container_type globals;
    std::pair<bool, MatchObject::context_type>
        retval(__rrmatch(rr, reactants, globals, i, j));
    // for (MatchObject::context_type::variable_container_type::const_iterator
    //     i(retval.second.globals.begin()); i != retval.second.globals.end(); ++i)
    // {
    //     std::cout << (*i).first << "=" << (*i).second << ";";
    // }
    // std::cout << std::endl;
    return retval.first;
}

// bool is_correspondent(const UnitSpecies& usp1, const UnitSpecies& usp2)
// {
//     if (usp1.name() != usp2.name() || usp1.num_sites() != usp2.num_sites())
//     {
//         return false;
//     }
// 
//     UnitSpecies::container_type::const_iterator i(usp1.begin()), j(usp2.begin());
//     while (i != usp1.end() && j != usp2.end())
//     {
//         if ((*i).first != (*j).first)
//         {
//             return false;
//         }
// 
//         ++i;
//         ++j;
//     }
//     return true;
// }

// void check_correspondences(const ReactionRule& rr)
// {
//     std::vector<unsigned int> correspondences;
//     unsigned int idx1(0), idx2(0);
// 
//     for (ReactionRule::reactant_container_type::const_iterator
//         i(rr.products().begin()); i != rr.products().end(); ++i)
//     {
//         for (Species::container_type::const_iterator
//             ii((*i).begin()); ii != (*i).end(); ++ii)
//         {
//             for (ReactionRule::reactant_container_type::const_iterator
//                 j(rr.reactants().begin()); j != rr.reactants().end(); ++j)
//             {
//                 for (Species::container_type::const_iterator
//                     jj((*j).begin()); jj != (*j).end(); ++jj)
//                 {
//                     if (is_correspondent(*ii, *jj))
//                     {
//                         ;
//                     }
//                     ++idx2;
//                 }
//             }
//             ++idx1;
//         }
//     }
// }

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
        std::pair<bool, MatchObject::context_type>
            retval(uspmatch(pttrn_, usp, ctx_));
        if (retval.first)
        {
            retval.second.iterators.push_back(pos);
            ++itr_;
            return retval;
        }
    }
    return std::make_pair(false, MatchObject::context_type());
}

} // ecell4
