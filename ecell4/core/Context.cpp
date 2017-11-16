#include "Context.hpp"
#include <string>
#include <sstream>


namespace ecell4
{

std::pair<bool, MatchObject::context_type> uspmatch(
    const UnitSpecies& pttrn, const UnitSpecies& usp,
    const MatchObject::context_type& org)
{
    std::pair<bool, MatchObject::context_type>
        retval(std::make_pair(false, org));
    MatchObject::context_type& ctx(retval.second);

    if (rbex::is_wildcard(pttrn.name()))
    {
        if (rbex::is_pass_wildcard(pttrn.name()))
        {
            throw NotSupported(
                "A pass wildcard '_0' is not allowed to be a name of Species.");
        }
        else if (rbex::is_named_wildcard(pttrn.name()))
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
                else if (rbex::is_pass_wildcard((*j).second.first))
                {
                    throw NotSupported(
                        "A pass wildcard '_0' is not allowed to be a state.");
                }
                else if (rbex::is_unnamed_wildcard((*j).second.first))
                {
                    ; // do nothing
                }
                else if (rbex::is_named_wildcard((*j).second.first))
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

            if (rbex::is_pass_wildcard((*j).second.second))
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
                else if (rbex::is_unnamed_wildcard((*j).second.second))
                {
                    continue;
                }
                else if (rbex::is_named_wildcard((*j).second.second))
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

    utils::get_mapper_mf<std::string, std::string>::type bond_cache;

    for (Species::container_type::const_iterator j(units2.begin());
        j != units2.end(); ++j)
    {
        units1.push_back(*j);

        for (UnitSpecies::container_type::const_iterator
            k((*j).begin()); k != (*j).end(); ++k)
        {
            const std::string& bond((*k).second.second);
            if (bond != "" && !rbex::is_wildcard(bond))
            {
                utils::get_mapper_mf<std::string, std::string>::type::const_iterator
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
        const unsigned int _gid = tag_units(groups, group_id, *i, connections, notyet);
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

std::vector<Species> ReactionRuleExpressionMatcher::generate()
{
    const unit_group_type res = this->genunits(this->compile());
    return group_units(res.units, res.groups, res.num_groups);
}

ReactionRuleExpressionMatcher::operation_type ReactionRuleExpressionMatcher::compile()
{
    typedef std::vector<UnitSpecies>::size_type size_type;
    typedef std::vector<UnitSpecies>::const_iterator const_iterator;

    operation_type res = {};
    std::vector<UnitSpecies>& products = res.products;
    std::vector<size_type>& correspo = res.correspo;
    std::vector<size_type>& removed = res.removed;

    // 1. Concatenate units of a pattern

    std::vector<UnitSpecies> reactants;
    for (ReactionRule::reactant_container_type::const_iterator
        i(pttrn_.reactants().begin()); i != pttrn_.reactants().end(); ++i)
    {
        std::vector<UnitSpecies> const units = (*i).units();
        reactants.reserve(reactants.size() + units.size());
        std::copy(units.begin(), units.end(), std::back_inserter(reactants));
    }

    res.reserved = reactants.size();

    int product_bond_stride = 0;
    for (ReactionRule::reactant_container_type::const_iterator
        i(pttrn_.products().begin()); i != pttrn_.products().end(); ++i)
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

ReactionRuleExpressionMatcher::unit_group_type ReactionRuleExpressionMatcher::genunits(const ReactionRuleExpressionMatcher::operation_type& operations)
{
    typedef std::vector<UnitSpecies>::size_type size_type;
    typedef std::vector<UnitSpecies>::const_iterator const_iterator;

    if (itr_ != matchers_.end())
    {
        // Failed to match
        return unit_group_type();;
    }

    const std::vector<UnitSpecies>& products = operations.products;
    const std::vector<size_type>& correspo = operations.correspo;
    const std::vector<size_type>& removed = operations.removed;
    const size_type& reserved = operations.reserved;

    unit_group_type res = {};
    std::vector<UnitSpecies>& units = res.units;
    std::vector<unsigned int>& groups = res.groups;
    unsigned int& num_groups = res.num_groups;

    // 3. Concatenate units given as reactants

    const context_type ctx(context());

    int bond_stride = 0;

    //XXX: for (std::vector<reactant_container_type::size_type>::const_iterator
    //XXX:     i(permutation_.begin()); i != permutation_.end(); ++i)
    //XXX: {
    //XXX:     bond_stride += concatenate_units(units, target_[*i], bond_stride);
    //XXX: }
    for (ReactionRule::reactant_container_type::const_iterator
        i(target_.begin()); i != target_.end(); ++i)
    {
        bond_stride += concatenate_units(units, *i, bond_stride);
    }

    // 4. Modify units

    utils::get_mapper_mf<std::string, std::string>::type bond_cache;

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

        for (std::vector<std::pair<size_type, size_type> >::const_iterator itr1(priorities.begin());
            itr1 != priorities.end(); ++itr1)
        {
            const UnitSpecies& op = products[(*itr1).second];
            const size_type& tgt = (*itr1).first;

            if (tgt >= units.size())
            {
                // 4-1. Create a new unit
                assert(tgt == units.size());
                // tgt = units.size();
                units.push_back(op);
                if (rbex::is_named_wildcard(op.name()))
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
                else if (rbex::is_wildcard((*i).second.first))
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
                else if (rbex::is_wildcard((*i).second.second))
                {
                    if (rbex::is_named_wildcard((*i).second.second))
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
                    utils::get_mapper_mf<std::string, std::string>::type::const_iterator
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

    utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type bondinfo;
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
                else if (rbex::is_wildcard(bond))
                {
                    std::stringstream ss;
                    ss << "A bond in a product is a wildcard ["
                        << (*i).name() << "(" << (*j).first << "^" << bond << ")].";
                    throw IllegalState(ss.str());
                }

                utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::iterator
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

    const ReactionRule::policy_type& policy(pttrn_.policy());

    if ((policy & ReactionRule::STRICT) && !(policy & (ReactionRule::DESTROY | ReactionRule::IMPLICIT)))
    {
        for (utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::const_iterator
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

    if (policy & ReactionRule::IMPLICIT)
    {
        for (std::vector<UnitSpecies>::iterator i(units.begin());
            i != units.end(); ++i)
        {
            for (UnitSpecies::container_type::size_type j(0);
                j != (*i).num_sites(); ++j)
            {
                UnitSpecies::container_type::value_type& site((*i).at(j));
                const std::string bond(site.second.second);
                if (bond != "" && !rbex::is_wildcard(bond) && bondinfo[bond].second != done)
                {
                    site.second.second = "";
                }
            }
        }
    }

    if (policy & ReactionRule::DESTROY)
    {
        for (utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::const_iterator
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
    const std::vector<unsigned int>& groups, const unsigned int num_groups)
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
        utils::get_mapper_mf<std::string, std::string>::type new_bonds;
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
                j != usp.num_sites(); ++j)
            {
                UnitSpecies::container_type::value_type&
                    site(usp.at(j));
                const std::string bond(site.second.second);
                if (bond == "" || rbex::is_wildcard(bond))
                {
                    continue;
                }

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

            products[idx].add_unit(usp);
        }

        // products[idx] = format_species(products[idx]);
    }

    products.erase(std::remove(products.begin(), products.end(), Species()), products.end());

    return products;
}

std::vector<ReactionRule> ReactionRuleExpressionMatcher::gen(const ReactionRule::reactant_container_type& reactants)
{
    typedef std::vector<ReactionRule> return_type;

    if (!this->match(reactants))
    {
        return return_type(0);
    }
    else if (pttrn_.reactants().size() == 0)
    {
        return return_type(
            1, ReactionRule(reactants, pttrn_.products(), pttrn_.k()));  // Zeroth-order reactions
    }

    std::vector<std::vector<UnitSpecies> > candidates;
    const operation_type op = this->compile();

    return_type res;

    do
    {
        const unit_group_type _res = this->genunits(op);

        std::vector<std::vector<UnitSpecies> >::iterator i(std::find(candidates.begin(), candidates.end(), _res.units));
        if (i != candidates.end())
        {
            ; // (*i).set_k((*i).k() + rr.k());
        }
        else
        {
            candidates.push_back(_res.units);
            res.push_back(ReactionRule(reactants, group_units(_res.units, _res.groups, _res.num_groups), pttrn_.k()));
        }
    }
    while (this->next());

    // return_type res;
    // res.reserve(candidates.size());
    // for (std::vector<std::vector<UnitSpecies> >::const_iterator i(candidates.begin());
    //     i != candidates.end(); ++i)
    // {
    //     res.push_back(ReactionRule(reactants, group_units(*i, pttrn_.policy()), pttrn_.k()));
    // }
    return res;
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

} // ecell4
