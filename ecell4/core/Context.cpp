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
    return (name.size() > 1 && name[0] == '_' && !is_pass_wildcard(name));
}

bool is_pass_wildcard(const std::string& name)
{
    return name == "_0";
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
        if (is_pass_wildcard(pttrn.name()))
        {
            throw NotSupported(
                "A pass wildcard '_0' is not allowed to be a name of Species.");
        }
        else if (is_named_wildcard(pttrn.name()))
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

            if (is_pass_wildcard((*j).second.second))
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

unsigned int __tag_units(
    std::vector<unsigned int>& retval,
    const unsigned int& group_id, const unsigned int& idx,
    const std::vector<UnitSpecies>& units,
    const std::vector<std::vector<std::vector<UnitSpecies>::size_type> >& adj)
{
    if (retval[idx] != units.size())
    {
        // assert(retval[idx] < group_id);
        return group_id;
    }

    retval[idx] = group_id;
    for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
        i(adj[idx].begin()); i != adj[idx].end(); ++i)
    {
        __tag_units(retval, group_id, *i, units, adj);
    }

    return group_id + 1;
}

std::pair<std::vector<unsigned int>, unsigned int> tag_units(
    const std::vector<UnitSpecies>& units,
    const std::vector<std::vector<std::vector<UnitSpecies>::size_type> >& adj)
{
    std::pair<std::vector<unsigned int>, unsigned int> retval;
    retval.first.resize(units.size(), units.size());
    retval.second = 0;
    for (unsigned int idx(0); idx != units.size(); ++idx)
    {
        // std::cout << idx << " connects with ";
        // for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
        //     i(adj[idx].begin()); i != adj[idx].end(); ++i)
        // {
        //     std::cout << *i;
        // }
        // std::cout << std::endl;

        retval.second = __tag_units(
            retval.first, retval.second, idx, units, adj);
    }

    return retval;
}

void check_correspondences(
    const std::vector<UnitSpecies>& reactant_units,
    const std::vector<UnitSpecies>& product_units,
    std::vector<std::vector<UnitSpecies>::size_type>& correspo,
    std::vector<std::vector<UnitSpecies>::size_type>& removed)
{
    correspo.clear();
    removed.clear();

    std::vector<UnitSpecies>::size_type num_units(reactant_units.size());
    std::vector<UnitSpecies>::size_type idx1(0), idx2(0);
    for (std::vector<UnitSpecies>::const_iterator i(product_units.begin());
        i != product_units.end(); ++i, ++idx1)
    {
        idx2 = 0;
        for (std::vector<UnitSpecies>::const_iterator j(reactant_units.begin());
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

    for (std::vector<UnitSpecies>::size_type i(0);
        i < reactant_units.size(); ++i)
    {
        if (std::find(correspo.begin(), correspo.end(), i)
            == correspo.end())
        {
            removed.push_back(i);
        }
    }
}

int concatenate_units(std::vector<UnitSpecies>& units1, const Species& sp, const int bond_stride)
{
    const std::vector<UnitSpecies> units2 = sp.units();
    units1.reserve(units1.size() + units2.size());

    int bond_ministride = 0;
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
                const int bondi = atoi(bond.c_str());
                // assert(bondi > 0);
                const std::string newbond = itos(bondi + bond_stride);
                units1.back().at(std::distance((*j).begin(), k)).second.second = newbond;
                bond_ministride = std::max(bond_ministride, bondi);
            }
        }
    }
    return bond_ministride;
}

std::vector<Species> group_units(
    const std::vector<UnitSpecies>& units, const ReactionRule::policy_type& policy)
{
    const unsigned int maxidx = units.size();
    utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type tmp;
    std::vector<std::vector<std::vector<UnitSpecies>::size_type> > adj;
    adj.resize(maxidx);

    for (std::vector<UnitSpecies>::const_iterator i(units.begin());
        i != units.end(); ++i)
    {
        const unsigned int idx(std::distance(units.begin(), i));

        for (UnitSpecies::container_type::const_iterator j((*i).begin());
            j != (*i).end(); ++j)
        {
            const std::string bond((*j).second.second);
            if (bond == "" || is_wildcard(bond))
            {
                continue;
            }

            utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::iterator
                itr(tmp.find(bond));
            if (itr == tmp.end())
            {
                tmp[bond] = std::make_pair((*j).first, idx);
            }
            else
            {
                if (tmp[bond].second == maxidx)
                {
                    ; //WARN: a duplicated bond found
                }

                adj[idx].push_back((*itr).second.second);
                adj[(*itr).second.second].push_back(idx);
                tmp[bond].second = units.size();  // This means the bond is already assigned.
            }
        }
    }

    if ((policy & ReactionRule::STRICT) && !(policy & (ReactionRule::DESTROY | ReactionRule::IMPLICIT)))
    {
        for (utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::const_iterator
            i(tmp.begin()); i != tmp.end(); ++i)
        {
            if ((*i).second.second != maxidx)
            {
                throw IllegalState("A bond is not resolved.");
            }
        }
    }

    std::pair<std::vector<unsigned int>, unsigned int>
        group_ids_pair(tag_units(units, adj));

    std::vector<Species> products;
    products.resize(group_ids_pair.second);
    // for (std::vector<UnitSpecies>::iterator i(units.begin());
    //     i != units.end(); ++i)
    // {
    //     products[group_ids_pair.first[std::distance(units.begin(), i)]].add_unit(*i);
    // }
    for (unsigned int idx(0); idx != group_ids_pair.second; ++idx)
    {
        utils::get_mapper_mf<std::string, std::string>::type new_bonds;
        unsigned int stride(1);

        for (std::vector<unsigned int>::iterator
            i(group_ids_pair.first.begin()); i != group_ids_pair.first.end(); ++i)
        {
            if (idx != *i)
            {
                continue;
            }

            UnitSpecies usp(
                units[std::distance(group_ids_pair.first.begin(), i)]); //XXX: copy
            for (UnitSpecies::container_type::size_type j(0);
                j != usp.num_sites(); ++j)
            {
                UnitSpecies::container_type::value_type&
                    site(usp.at(j));
                const std::string bond(site.second.second);
                if (bond == "" || is_wildcard(bond))
                {
                    continue;
                }

                if ((policy & ReactionRule::IMPLICIT) && tmp[bond].second != maxidx)
                {
                    site.second.second = "";
                }
                else
                {
                    utils::get_mapper_mf<std::string, std::string>::type::const_iterator
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
            }

            products[idx].add_unit(usp);
        }

        // products[idx] = format_species(products[idx]);
    }

    if (policy & ReactionRule::DESTROY)
    {
        std::vector<unsigned int> removed;
        for (utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::const_iterator
            i(tmp.begin()); i != tmp.end(); ++i)
        {
            if ((*i).second.second != maxidx)
            {
                removed.push_back(group_ids_pair.first[(*i).second.second]);
            }
        }
        std::sort(removed.begin(), removed.end());
        removed.erase(std::unique(removed.begin(), removed.end()), removed.end());

        for (std::vector<unsigned int>::const_reverse_iterator
            i(removed.rbegin()); i != removed.rend(); ++i)
        {
            products.erase(products.begin() + *i);
        }
    }

    return products;
}

std::vector<Species> ReactionRuleExpressionMatcher::generate()
{
    if (itr_ != matchers_.end())
    {
        return std::vector<Species>();
    }
    else if (pttrn_.reactants().size() == 0)
    {
        return pttrn_.products();  // XXX: zero-th order reaction
    }

    const context_type ctx(context());

    std::vector<UnitSpecies> reactant_units, product_units;

    for (ReactionRule::reactant_container_type::const_iterator
        i(pttrn_.reactants().begin()); i != pttrn_.reactants().end(); ++i)
    {
        std::vector<UnitSpecies> units = (*i).units();
        reactant_units.reserve(reactant_units.size() + units.size());
        std::copy(units.begin(), units.end(), std::back_inserter(reactant_units));
    }

    for (ReactionRule::reactant_container_type::const_iterator
        i(pttrn_.products().begin()); i != pttrn_.products().end(); ++i)
    {
        std::vector<UnitSpecies> units = (*i).units();
        product_units.reserve(product_units.size() + units.size());
        std::copy(units.begin(), units.end(), std::back_inserter(product_units));
    }

    std::vector<std::vector<UnitSpecies>::size_type> correspo, removed;
    check_correspondences(reactant_units, product_units, correspo, removed);

    int bond_stride = 0;
    std::vector<UnitSpecies> units;
    // for (ReactionRule::reactant_container_type::const_iterator
    //     i(target_.begin()); i != target_.end(); ++i)
    // {
    //     bond_stride += concatenate_units(units, *i, bond_stride);
    // }
    for (std::vector<reactant_container_type::size_type>::const_iterator
        i(permutation_.begin()); i != permutation_.end(); ++i)
    {
        bond_stride += concatenate_units(units, target_[*i], bond_stride);
    }
    // for (context_type::iterator_container_type::const_iterator
    //     i(ctx.iterators.begin()); i != ctx.iterators.end(); ++i)
    // {
    //     std::cout << (*i);
    // }
    // std::cout << std::endl;

    std::vector<UnitSpecies>::size_type idx1(0);
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
            if (is_named_wildcard((*itr1).name()))
            {
                context_type::variable_container_type::const_iterator
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
            else if (is_wildcard((*i).second.first))
            {
                if ((*i).second.first.size() != 1)
                {
                    context_type::variable_container_type::const_iterator
                        itr(ctx.globals.find((*i).second.first));
                    if (itr == ctx.globals.end())
                    {
                        ; //XXX: an invalid global name given
                        std::cerr << "An invalid global name given ["
                            << (*i).second.first << "]" << std::endl;
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

            if ((*i).second.second == "")
            {
                site.second.second = "";
            }
            else if (is_wildcard((*i).second.second))
            {
                if ((*i).second.second.size() != 1)
                {
                    ; //XXX: no named wildcard is allowed here
                    ; //XXX: how about the pass wildcard?
                }

                ; // do nothing
            }
            else
            {
                unsigned int stride(0), label(0);
                for (ReactionRule::product_container_type::const_iterator
                    j(pttrn_.products().begin()); j != pttrn_.products().end(); ++j)
                {
                    stride += (*j).units().size();
                    if (stride > idx1)
                    {
                        label = atoi((*i).second.second.c_str()) * pttrn_.products().size()
                            + std::distance(pttrn_.products().begin(), j);
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
                    ++bond_stride;
                    site.second.second = itos(bond_stride);
                    bond_cache[label] = site.second.second;
                }
            }

            ++i;
        }
    }

    std::vector<std::vector<UnitSpecies>::size_type> removed_new;
    for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
        i(removed.begin()); i != removed.end(); ++i)
    {
        removed_new.push_back(ctx.iterators[(*i)]);
    }
    std::sort(removed_new.begin(), removed_new.end());

    for (std::vector<std::vector<UnitSpecies>::size_type>::const_reverse_iterator
        i(removed_new.rbegin()); i != removed_new.rend(); ++i)
    {
        units.erase(units.begin() + *i);
    }

    return group_units(units, pttrn_.policy());
    // std::cout << std::endl << "before: ";
    // for (ReactionRule::reactant_container_type::const_iterator i(target_.begin());
    //     i != target_.end(); ++i)
    // {
    //     std::cout << (*i).serial() << ".";
    // }
    // std::cout << std::endl;
    // std::cout << "after: ";
    // for (std::vector<UnitSpecies>::const_iterator i(units.begin());
    //     i != units.end(); ++i)
    // {
    //     std::cout << (*i).serial() << ".";
    // }
    // std::cout << std::endl;
    // std::cout << "after: ";
    // for (std::vector<unsigned int>::const_iterator i(group_ids_pair.first.begin());
    //     i != group_ids_pair.first.end(); ++i)
    // {
    //     std::cout << (*i) << ".";
    // }
    // std::cout << std::endl;
    // std::cout << "products: ";
    // for (std::vector<Species>::const_iterator i(products.begin());
    //     i != products.end(); ++i)
    // {
    //     std::cout << (*i).serial() << ", ";
    // }
    // std::cout << std::endl;
    // return products;
}

std::pair<bool, MatchObject::context_type> MatchObject::next()
{
    std::vector<UnitSpecies>::const_iterator itr_start = target_.begin();
    for (; itr_ != target_.end(); ++itr_)
    {
        const Species::container_type::difference_type
            pos(distance(itr_start, itr_));
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

// ReactionRule create_reaction_rule_formatted(
//     const ReactionRule::reactant_container_type& reactants,
//     const ReactionRule::product_container_type& products, const Real k)
// {
//     ReactionRule rr(reactants, ReactionRule::product_container_type(), k);
//     for (ReactionRule::product_container_type::const_iterator i(products.begin());
//         i != products.end(); ++i)
//     {
//         rr.add_product(format_species(*i));
//     }
//     return rr;
// }

} // ecell4
