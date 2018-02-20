#include <algorithm>

#include "exceptions.hpp"
#include "NetfreeModel.hpp"


namespace ecell4
{

std::vector<ReactionRule> NetfreeModel::query_reaction_rules(
    const Species& sp) const
{
    ReactionRule::reactant_container_type reactants(1, sp);
    std::vector<ReactionRule> retval;
    for (reaction_rule_container_type::const_iterator i(reaction_rules_.begin());
        i != reaction_rules_.end(); ++i)
    {
        const std::vector<ReactionRule> generated = (*i).generate(reactants);
        // retval.insert(retval.end(), generated.begin(), generated.end());
        retval.reserve(retval.size() + generated.size());
        for (std::vector<ReactionRule>::const_iterator j(generated.begin());
            j != generated.end(); ++j)
        {
            // const ReactionRule rr = create_reaction_rule_formatted(*j);
            const ReactionRule rr = format_reaction_rule_with_nosort(*j);
            std::vector<ReactionRule>::iterator
                it = std::find(retval.begin(), retval.end(), rr);
            if (it == retval.end())
            {
                retval.push_back(rr);
            }
            else
            {
                (*it).set_k((*it).k() + rr.k());
            }
        }
    }
    return retval;
}

struct reaction_rule_product_unary_predicator
{
    typedef ReactionRule element_type;

    reaction_rule_product_unary_predicator(const element_type& target)
        : target_(target)
    {
        ; // do nothing
    }

    bool operator()(const element_type& v)
    {
        return v.products() == target_.products();
    }

protected:

    element_type target_;
};

std::vector<ReactionRule> generate_reaction_rules(
    const ReactionRule& org, const Species& sp1, const Species& sp2)
{
    if (org.reactants().size() != 2)
    {
        return std::vector<ReactionRule>(0);
    }

    ReactionRule::reactant_container_type reactants(2);
    reactants[0] = sp1;
    reactants[1] = sp2;
    std::vector<ReactionRule> res = org.generate(reactants);

    if (org.reactants()[0] != org.reactants()[1])
    {
        reactants[0] = sp2;
        reactants[1] = sp1;
        std::vector<ReactionRule> _res = org.generate(reactants);
        std::copy(_res.begin(), _res.end(), back_inserter(res));
    }
    return res;
}

std::vector<ReactionRule> NetfreeModel::query_reaction_rules(
    const Species& sp1, const Species& sp2) const
{
    std::vector<ReactionRule> retval;
    for (reaction_rule_container_type::const_iterator i(reaction_rules_.begin());
        i != reaction_rules_.end(); ++i)
    {
        const std::vector<ReactionRule> generated = generate_reaction_rules(*i, sp1, sp2);
        // retval.insert(retval.end(), generated.begin(), generated.end());
        retval.reserve(retval.size() + generated.size());
        for (std::vector<ReactionRule>::const_iterator j(generated.begin());
            j != generated.end(); ++j)
        {
            // const ReactionRule rr = create_reaction_rule_formatted(*j);
            const ReactionRule rr = format_reaction_rule_with_nosort(*j);
            std::vector<ReactionRule>::iterator
                it = std::find(retval.begin(), retval.end(), rr);
            if (it == retval.end())
            {
                retval.push_back(rr);
            }
            else
            {
                (*it).set_k((*it).k() + rr.k());
            }
        }
    }

    if (effective_)
    {
        for (std::vector<ReactionRule>::iterator i(retval.begin()); i != retval.end(); ++i)
        {
            const ReactionRule& rr(*i);
            if (rr.reactants()[0] == rr.reactants()[1])
            {
                (*i).set_k(rr.k() * 0.5);
            }
        }
    }
    return retval;
}

Integer NetfreeModel::apply(const Species& pttrn, const Species& sp) const
{
    return pttrn.count(sp);
}

std::vector<ReactionRule> NetfreeModel::apply(
    const ReactionRule& rr, const ReactionRule::reactant_container_type& reactants) const
{
    return rr.generate(reactants);
}

void NetfreeModel::add_species_attribute(const Species& sp)
{
    if (has_species_attribute_exact(sp))
    {
        throw AlreadyExists("species already exists");
    }
    species_attributes_.push_back(sp);
}

void NetfreeModel::remove_species_attribute(const Species& sp)
{
    species_container_type::iterator i(
        std::find(species_attributes_.begin(), species_attributes_.end(), sp));
    if (i == species_attributes_.end())
    {
        std::ostringstream message;
        message << "Speices [" << sp.serial() << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }
    species_attributes_.erase(i);
}

bool NetfreeModel::has_species_attribute(const Species& sp) const
{
    return has_species_attribute_exact(sp);
}

bool NetfreeModel::has_species_attribute_exact(const Species& sp) const
{
    species_container_type::const_iterator i(
        std::find(species_attributes_.begin(), species_attributes_.end(), sp));
    return (i != species_attributes_.end());
}

void NetfreeModel::add_reaction_rule(const ReactionRule& rr)
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    if (i != reaction_rules_.end())
    {
        throw AlreadyExists("reaction rule already exists");
    }

    // const reaction_rule_container_type::size_type idx(reaction_rules_.size());
    reaction_rules_.push_back(rr);

    // if (rr.reactants().size() == 1)
    // {
    //     first_order_reaction_rules_map_[rr.reactants()[0].serial()].push_back(idx);
    // }
    // else if (rr.reactants().size() == 2)
    // {
    //     const Species::serial_type
    //         serial1(rr.reactants()[0].serial()),
    //         serial2(rr.reactants()[1].serial());
    //     const std::pair<Species::serial_type, Species::serial_type>
    //         key(serial1 < serial2?
    //             std::make_pair(serial1, serial2):
    //             std::make_pair(serial2, serial1));
    //     second_order_reaction_rules_map_[key].push_back(idx);
    // }
    // else
    // {
    //     ;
    // }
}

void NetfreeModel::remove_reaction_rule(const ReactionRule& rr)
{
    reaction_rule_container_type::iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    if (i == reaction_rules_.end())
    {
        throw NotFound("reaction rule not found");
    }
    reaction_rules_.erase(i);

    // reaction_rule_container_type::size_type const
    //     idx(i - reaction_rules_.begin()), last_idx(reaction_rules_.size() - 1);
    // if (rr.reactants().size() == 1)
    // {
    //     first_order_reaction_rules_map_type::iterator
    //         j(first_order_reaction_rules_map_.find(rr.reactants()[0].serial()));
    //     if (j == first_order_reaction_rules_map_.end())
    //     {
    //         throw IllegalState("no corresponding map key found");
    //     }

    //     first_order_reaction_rules_map_type::mapped_type::iterator
    //         k(std::remove((*j).second.begin(), (*j).second.end(), idx));
    //     if (k == (*j).second.end())
    //     {
    //         throw IllegalState("no corresponding map value found");
    //     }
    //     else
    //     {
    //         (*j).second.erase(k, (*j).second.end());
    //     }
    // }
    // else if (rr.reactants().size() == 2)
    // {
    //     second_order_reaction_rules_map_type::iterator
    //         j(second_order_reaction_rules_map_.find(std::make_pair(
    //             rr.reactants()[0].serial(), rr.reactants()[1].serial())));
    //     if (j == second_order_reaction_rules_map_.end())
    //     {
    //         throw IllegalState("no corresponding map key found");
    //     }

    //     second_order_reaction_rules_map_type::mapped_type::iterator
    //         k(std::remove((*j).second.begin(), (*j).second.end(), idx));
    //     if (k == (*j).second.end())
    //     {
    //         throw IllegalState("no corresponding map value found");
    //     }
    //     else
    //     {
    //         (*j).second.erase(k, (*j).second.end());
    //     }
    // }

    // if (idx < last_idx)
    // {
    //     reaction_rule_container_type::value_type const
    //         last_value(reaction_rules_[last_idx]);
    //     (*i) = last_value;

    //     if (last_value.reactants().size() == 1)
    //     {
    //         first_order_reaction_rules_map_type::iterator
    //             j(first_order_reaction_rules_map_.find(
    //                 last_value.reactants()[0].serial()));
    //         if (j == first_order_reaction_rules_map_.end())
    //         {
    //             throw IllegalState("no corresponding map key for the last found");
    //         }

    //         first_order_reaction_rules_map_type::mapped_type::iterator
    //             k(std::remove((*j).second.begin(), (*j).second.end(), last_idx));
    //         if (k == (*j).second.end())
    //         {
    //             throw IllegalState("no corresponding map value found");
    //         }
    //         else
    //         {
    //             (*j).second.erase(k, (*j).second.end());
    //         }
    //         (*j).second.push_back(idx);
    //     }
    //     else if (last_value.reactants().size() == 2)
    //     {
    //         second_order_reaction_rules_map_type::iterator
    //             j(second_order_reaction_rules_map_.find(std::make_pair(
    //                 last_value.reactants()[0].serial(),
    //                 last_value.reactants()[1].serial())));
    //         if (j == second_order_reaction_rules_map_.end())
    //         {
    //             throw IllegalState("no corresponding map key for the last found");
    //         }
    //         second_order_reaction_rules_map_type::mapped_type::iterator
    //             k(std::remove((*j).second.begin(), (*j).second.end(), last_idx));
    //         if (k == (*j).second.end())
    //         {
    //             throw IllegalState("no corresponding map value found");
    //         }
    //         else
    //         {
    //             (*j).second.erase(k, (*j).second.end());
    //         }
    //         (*j).second.push_back(idx);
    //     }
    // }

    // reaction_rules_.pop_back();
}

bool NetfreeModel::has_reaction_rule(const ReactionRule& rr) const
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    return (i != reaction_rules_.end());
}

boost::shared_ptr<Model> NetfreeModel::expand(
    const std::vector<Species>& sp, const Integer max_itr,
    const std::map<Species, Integer>& max_stoich) const
{
    return extras::generate_network_from_netfree_model(
        *this, sp, max_itr, max_stoich).first;
}


boost::shared_ptr<Model> NetfreeModel::expand(
    const std::vector<Species>& sp, const Integer max_itr) const
{
    return extras::generate_network_from_netfree_model(
        *this, sp, max_itr).first;
}

boost::shared_ptr<Model> NetfreeModel::expand(
    const std::vector<Species>& sp) const
{
    const Integer max_itr(30);
    std::pair<boost::shared_ptr<NetworkModel>, bool>
        retval(extras::generate_network_from_netfree_model(*this, sp, max_itr));
    if (retval.second)
    {
        return retval.first;
    }
    else
    {
        return boost::shared_ptr<NetworkModel>(); // return null
    }
}

namespace extras
{

bool check_stoichiometry(const Species& sp,
    const std::map<Species, Integer>& max_stoich)
{
    for (std::map<Species, Integer>::const_iterator i(max_stoich.begin());
        i != max_stoich.end(); ++i)
    {
        if ((*i).first.count(sp) > (*i).second)
        {
            return false;
        }
    }
    return true;
}

bool check_stoichiometry(const ReactionRule& rr,
    const std::map<Species, Integer>& max_stoich)
{
    for (ReactionRule::product_container_type::const_iterator
        i(rr.products().begin()); i != rr.products().end(); ++i)
    {
        if (!check_stoichiometry(*i, max_stoich))
        {
            return false;
        }
    }
    return true;
}

void __add_reaction_rules(
    const std::vector<ReactionRule>& reaction_rules,
    std::vector<ReactionRule>& reactions, std::vector<Species>& newseeds,
    const std::vector<Species>& seeds,
    const std::map<Species, Integer>& max_stoich)
{
    for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
        i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        if (!check_stoichiometry(rr, max_stoich))
        {
            continue;
        }

        reactions.push_back(rr);

        for (ReactionRule::product_container_type::const_iterator
            j(rr.products().begin()); j != rr.products().end(); ++j)
        {
            const Species sp(format_species(*j));
            if (std::find(newseeds.begin(), newseeds.end(), sp)
                == newseeds.end()
                && std::find(seeds.begin(), seeds.end(), sp)
                == seeds.end())
            {
                newseeds.push_back(sp);
            }
        }
    }
}

void __generate_recurse(
    const NetfreeModel& nfm, std::vector<ReactionRule>& reactions,
    std::vector<Species>& seeds1, std::vector<Species>& seeds2,
    const std::map<Species, Integer>& max_stoich)
{
    std::vector<Species> newseeds;
    seeds2.insert(seeds2.begin(), seeds1.begin(), seeds1.end());

    for (NetfreeModel::reaction_rule_container_type::const_iterator
        i(nfm.reaction_rules().begin()); i != nfm.reaction_rules().end(); ++i)
    {
        const ReactionRule& rr(*i);

        switch (rr.reactants().size())
        {
        case 0:
            continue;
        case 1:
            for (std::vector<Species>::const_iterator j(seeds1.begin());
                j != seeds1.end(); ++j)
            {
                ReactionRule::reactant_container_type reactants(1);
                reactants[0] = *j;
                __add_reaction_rules(
                    rr.generate(reactants), reactions, newseeds, seeds2, max_stoich);
            }
            break;
        case 2:
            for (std::vector<Species>::const_iterator j(seeds1.begin());
                j != seeds1.end(); ++j)
            {
                const std::vector<Species>::const_iterator start(
                    seeds2.begin()
                    + std::distance<std::vector<Species>::const_iterator>(
                        seeds1.begin(), j));
                for (std::vector<Species>::const_iterator
                    k(start); k != seeds2.end(); ++k)
                {
                    __add_reaction_rules(
                        generate_reaction_rules(rr, *j, *k),
                        reactions, newseeds, seeds2, max_stoich);
                }
            }
            break;
        default:
            throw NotImplemented(
                "No reaction rule with more than two reactants is accepted.");
        }
    }

    seeds1.swap(newseeds);
}

std::pair<boost::shared_ptr<NetworkModel>, bool> generate_network_from_netfree_model(
    const NetfreeModel& nfm, const std::vector<Species>& seeds, const Integer max_itr,
    const std::map<Species, Integer>& max_stoich)
{
    std::vector<ReactionRule> reactions;
    std::vector<Species> seeds1(seeds);
    std::vector<Species> seeds2;

    for (NetfreeModel::reaction_rule_container_type::const_iterator
        i(nfm.reaction_rules().begin()); i != nfm.reaction_rules().end(); ++i)
    {
        const ReactionRule& rr(*i);
        if (rr.reactants().size() == 0 && check_stoichiometry(rr, max_stoich))
        {
            reactions.push_back(rr);
            for (ReactionRule::product_container_type::const_iterator
                j(rr.products().begin()); j != rr.products().end(); ++j)
            {
                const Species sp(format_species(*j));
                if (std::find(seeds1.begin(), seeds1.end(), sp)
                    == seeds1.end())
                {
                    seeds1.push_back(sp);
                }
            }
        }
    }

    Integer cnt(0);
    while (seeds1.size() > 0 && cnt < max_itr)
    {
        __generate_recurse(nfm, reactions, seeds1, seeds2, max_stoich);
        cnt += 1;
    }

    bool is_completed;
    if (seeds1.size() != 0)
    {
        is_completed = false;
        seeds2.insert(seeds2.begin(), seeds1.begin(), seeds1.end());
    }
    else
    {
        is_completed = true;
    }

    boost::shared_ptr<NetworkModel> nwm(new NetworkModel());
    for (std::vector<Species>::const_iterator i(seeds2.begin());
        i != seeds2.end(); ++i)
    {
        (*nwm).add_species_attribute(nfm.apply_species_attributes(*i));
    }
    if (nfm.effective())
    {
        for (std::vector<ReactionRule>::const_iterator i(reactions.begin());
            i != reactions.end(); ++i)
        {
            ReactionRule rr(format_reaction_rule(*i));
            if (rr.reactants().size() == 2 && rr.reactants()[0] == rr.reactants()[1])
            {
                rr.set_k(rr.k() * 0.5);
            }

            (*nwm).add_reaction_rule(rr);
        }
    }
    else
    {
        for (std::vector<ReactionRule>::const_iterator i(reactions.begin());
            i != reactions.end(); ++i)
        {
            (*nwm).add_reaction_rule(format_reaction_rule(*i));
        }
    }
    return std::make_pair(nwm, is_completed);
}

} // extras

} // ecell4
