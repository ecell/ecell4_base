#include "Context.hpp"
#include <string>


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

    ReactionRuleExpressionMatcher rrexp(rr);
    return rrexp.match(reactants);
    // if (rr.reactants().size() != reactants.size())
    // {
    //     return false;
    // }

    // ReactionRule::reactant_container_type::const_iterator
    //     i(rr.reactants().begin()), j(reactants.begin());
    // MatchObject::context_type::variable_container_type globals;
    // std::pair<bool, MatchObject::context_type>
    //     retval(__rrmatch(rr, reactants, globals, i, j));
    // // for (MatchObject::context_type::variable_container_type::const_iterator
    // //     i(retval.second.globals.begin()); i != retval.second.globals.end(); ++i)
    // // {
    // //     std::cout << (*i).first << "=" << (*i).second << ";";
    // // }
    // // std::cout << std::endl;
    // return retval.first;
}

Integer count_rrmatches(const ReactionRule& rr,
    const ReactionRule::reactant_container_type& reactants)
{
    ReactionRuleExpressionMatcher rrexp(rr);
    if (!rrexp.match(reactants))
    {
        return 0;
    }
    Integer n(1);
    while (rrexp.next())
    {
        ++n;
    }
    return n;
}

bool is_correspondent(const UnitSpecies& usp1, const UnitSpecies& usp2)
{
    if (usp1.name() != usp2.name() || usp1.num_sites() != usp2.num_sites())
    {
        return false;
    }

    UnitSpecies::container_type::const_iterator i(usp1.begin()), j(usp2.begin());
    while (i != usp1.end() && j != usp2.end())
    {
        if ((*i).first != (*j).first)
        {
            return false;
        }

        ++i;
        ++j;
    }
    return true;
}

std::string itos(unsigned int val)
{
    std::stringstream ss;
    ss << val;
    return ss.str();
}

void rrgenerate(const ReactionRule& rr,
    const ReactionRule::reactant_container_type& reactants)
{
    ReactionRuleExpressionMatcher rrexp(rr);
    if (!rrexp.match(reactants))
    {
        return;
    }
    const ReactionRuleExpressionMatcher::context_type ctx(rrexp.context());

    std::vector<UnitSpecies> reactant_units, product_units;

    for (ReactionRule::reactant_container_type::const_iterator
        i(rr.reactants().begin()); i != rr.reactants().end(); ++i)
    {
        reactant_units.reserve(reactant_units.size() + (*i).num_units());
        std::copy((*i).begin(), (*i).end(), std::back_inserter(reactant_units));
    }

    for (ReactionRule::reactant_container_type::const_iterator
        i(rr.products().begin()); i != rr.products().end(); ++i)
    {
        product_units.reserve(product_units.size() + (*i).num_units());
        std::copy((*i).begin(), (*i).end(), std::back_inserter(product_units));
    }

    std::vector<UnitSpecies>::size_type num_units(reactant_units.size());
    std::vector<std::vector<UnitSpecies>::size_type> correspo;
    std::vector<UnitSpecies>::size_type idx1(0), idx2(0);
    for (std::vector<UnitSpecies>::iterator i(product_units.begin());
        i != product_units.end(); ++i, ++idx1)
    {
        idx2 = 0;
        for (std::vector<UnitSpecies>::iterator j(reactant_units.begin());
            j != reactant_units.end(); ++j, ++idx2)
        {
            if (is_correspondent(*i, *j))
            {
                if (correspo.size() > idx1)
                {
                    ; //WARN: multiple correspondence found
                }
                else if (std::find(correspo.begin(), correspo.end(), idx2)
                    != correspo.end())
                {
                    ; //WARN: multiple correspondence skipped
                }
                else
                {
                    correspo.push_back(idx2);
                }
            }
        }

        if (correspo.size() == idx1)
        {
            correspo.push_back(num_units);
            ++num_units;
        }
    }

    std::vector<std::vector<UnitSpecies>::size_type> removed;
    for (std::vector<UnitSpecies>::size_type i(0);
        i < reactant_units.size(); ++i)
    {
        if (std::find(correspo.begin(), correspo.end(), i)
            == correspo.end())
        {
            removed.push_back(i);
        }
    }

    int bond_stride1(0), bond_stride2(0);
    std::vector<UnitSpecies> units;
    for (ReactionRule::reactant_container_type::const_iterator
        i(reactants.begin()); i != reactants.end(); ++i)
    {
        units.reserve(units.size() + (*i).num_units());
        // std::copy((*i).begin(), (*i).end(), std::back_inserter(units));

        for (Species::container_type::const_iterator j((*i).begin());
            j != (*i).end(); ++j)
        {
            units.push_back(*j);
            for (UnitSpecies::container_type::const_iterator
                k((*j).begin()); k != (*j).end(); ++k)
            {
                if ((*k).second.second != "" && (*k).second.second[0] != '_')
                {
                    units.back().at(std::distance((*j).begin(), k)).second.second = itos(atoi((*k).second.second.c_str()) + bond_stride1);
                    bond_stride2 = std::max(bond_stride2, atoi((*k).second.second.c_str()));
                }
            }
        }

        bond_stride1 += bond_stride2;
        bond_stride2 = 0;
    }

    // for (MatchObject::context_type::iterator_container_type::const_iterator
    //     i(ctx.iterators.begin()); i != ctx.iterators.end(); ++i)
    // {
    //     std::cout << (*i);
    // }
    // std::cout << std::endl;

    idx1 = 0;
    utils::get_mapper_mf<unsigned int, std::string>::type bond_cache;
    for (std::vector<UnitSpecies>::iterator itr1(product_units.begin());
        itr1 != product_units.end(); ++itr1, ++idx1)
    {
        std::vector<UnitSpecies>::size_type tgt(0);
        std::vector<UnitSpecies>::size_type idx2(correspo[idx1]);
        if (idx2 >= reactant_units.size())
        {
            tgt = units.size();
            units.push_back(*itr1);
            if ((*itr1).name().size() > 1 && (*itr1).name()[0] == '_')
            {
                MatchObject::context_type::variable_container_type::const_iterator
                    itr(ctx.globals.find((*itr1).name()));
                if (itr == ctx.globals.end())
                {
                    ; //XXX: an invalid global name given
                }
                else
                {
                    units.back().set_name((*itr).second);
                }
            }

            // continue;
        }
        else
        {
            tgt = ctx.iterators[idx2];
        }
        // std::cout << "(" << idx2 << "," << tgt << ")";

        UnitSpecies::container_type::const_iterator i((*itr1).begin());
        while (i != (*itr1).end())
        {
            UnitSpecies::container_type::value_type&
                site(units[tgt].at((*i).first));

            if ((*i).second.first == "")
            {
                ; // do nothing
            }
            else if ((*i).second.first[0] == '_')
            {
                if ((*i).second.first.size() != 1)
                {
                    MatchObject::context_type::variable_container_type::const_iterator
                        itr(ctx.globals.find((*i).second.first));
                    if (itr == ctx.globals.end())
                    {
                        ; //XXX: an invalid global name given
                    }
                    else
                    {
                        site.first = (*itr).second;
                    }
                }
            }
            else
            {
                site.second.first = (*i).second.first;
            }

            if ((*i).second.second == "")
            {
                site.second.second = "";
            }
            else if ((*i).second.second[0] == '_')
            {
                if ((*i).second.second.size() != 1)
                {
                    ; //XXX: no named wildcard is allowed here
                }

                ; // do nothing
            }
            else
            {
                unsigned int stride(0), label(0);
                for (ReactionRule::product_container_type::const_iterator
                    j(rr.products().begin()); j != rr.products().end(); ++j)
                {
                    stride += (*j).num_units();
                    if (stride > idx1)
                    {
                        label = atoi((*i).second.second.c_str()) * rr.products().size() + std::distance(rr.products().begin(), j);
                        break;
                    }
                }

                utils::get_mapper_mf<unsigned int, std::string>::type::const_iterator
                    itr(bond_cache.find(label));
                if (itr != bond_cache.end())
                {
                    site.second.second = (*itr).second;
                }
                else
                {
                    ++bond_stride1;
                    site.second.second = itos(bond_stride1);
                    bond_cache[label] = site.second.second;
                }
            }

            ++i;
        }
    }

    std::vector<std::vector<UnitSpecies>::size_type> removed_new;
    for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator i(removed.begin()); i != removed.end(); ++i)
    {
        removed_new.push_back(ctx.iterators[(*i)]);
    }
    std::sort(removed_new.begin(), removed_new.end());
    for (std::vector<std::vector<UnitSpecies>::size_type>::const_reverse_iterator i(removed_new.rbegin()); i != removed_new.rend(); ++i)
    {
        units.erase(units.begin() + *i);
    }

    //XXX: under construction

    std::cout << std::endl << "before: ";
    for (ReactionRule::reactant_container_type::const_iterator i(reactants.begin());
        i != reactants.end(); ++i)
    {
        std::cout << (*i).serial() << ".";
    }
    std::cout << std::endl;
    std::cout << "after:  ";
    for (std::vector<UnitSpecies>::const_iterator i(units.begin());
        i != units.end(); ++i)
    {
        std::cout << (*i).serial() << ".";
    }
    std::cout << std::endl;
}

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
