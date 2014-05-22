#include "Species.hpp"
#include "Context.hpp"

#include <algorithm>


namespace ecell4
{

bool Species::operator==(const Species& rhs) const
{
    return (serial() == rhs.serial());
}

bool Species::operator<(const Species& rhs) const
{
    return (serial() < rhs.serial());
}

bool Species::operator>(const Species& rhs) const
{
    return (serial() > rhs.serial());
}

bool Species::match(const Species& target) const
{
    return spmatch(*this, target);
}

std::string serialize_unit_species_masked(const UnitSpecies& usp)
{
    if (usp.num_sites() == 0)
    {
        return usp.name();
    }

    std::vector<std::string> unstated, stated;
    for (UnitSpecies::container_type::const_iterator i(usp.begin());
        i != usp.end(); ++i)
    {
        const std::string&
            state((*i).second.first), bond((*i).second.second);
        if (state.size() > 0)
        {
            stated.push_back(
                (*i).first + "="
                + (bond.size() > 0? state + "^_" : state));
        }
        else
        {
            unstated.push_back(
                bond.size() > 0? (*i).first + "^_" : (*i).first);
        }
    }

    std::sort(unstated.begin(), unstated.end());
    std::sort(stated.begin(), stated.end());
    return usp.name() + "(" + boost::algorithm::join(unstated, ",")
        + (unstated.size() > 0 && stated.size() > 0? "," : "")
        + boost::algorithm::join(stated, ",") + ")";
}

class unit_species_comparerator
{
public:

    // typedef Species::container_type::size_type index_type;
    typedef unsigned int index_type;
    typedef std::pair<index_type, std::string> site_type;
    typedef utils::get_mapper_mf<std::string, std::vector<site_type> >::type connection_container_type;

public:

    unit_species_comparerator(const Species& sp)
        : root_(sp)
    {
        initialize();
    }

    void initialize()
    {
        connections_.clear();
        for (index_type idx(0); idx < root_.num_units(); ++idx)
        {
            const UnitSpecies usp(root_.at(idx));
            for (UnitSpecies::container_type::const_iterator i(usp.begin());
                 i != usp.end(); ++i)
            {
                if ((*i).second.second == "" || (*i).second.second[0] == '_')
                {
                    continue;
                }

                if (connections_.find((*i).second.second) == connections_.end())
                {
                    connections_.insert(std::make_pair((*i).second.second, std::vector<site_type>()));
                }
                connections_[(*i).second.second].push_back(std::make_pair(idx, (*i).first));
            }
        }
    }

    int compare(const index_type& val1, const index_type& val2)
    {
        if (val1 == val2)
        {
            return 0;
        }

        const std::pair<index_type, index_type> pair_key(
            (val1 < val2)? std::make_pair(val1, val2) : std::make_pair(val1, val2));
        if (std::binary_search(ignores_.begin(), ignores_.end(), pair_key))
        {
            return 0;
        }

        const UnitSpecies& lhs(root_.units().at(val1)), rhs(root_.units().at(val2));

        if (lhs.name() != rhs.name())
        {
            return (lhs.name() < rhs.name()? 1 : -1);
        }

        UnitSpecies::container_type::const_iterator i(lhs.begin()), j(rhs.begin());
        while (i != lhs.end() && j != rhs.end())
        {
            if ((*i).first != (*j).first)
            {
                // std::cout << "[1] " << lhs.serial() << "(" << val1 << ") vs " << rhs.serial() << "(" << val2 << ") -> " << (*i).first << " < " << (*j).first << std::endl;
                return ((*i).first < (*j).first? 1 : -1);
            }
            else if ((*i).second.first != (*j).second.first)
            {
                // std::cout << "[2] " << lhs.serial() << "(" << val1 << ") vs " << rhs.serial() << "(" << val2 << ")" << std::endl;
                return ((*i).second.first < (*j).second.first? 1 : -1);
            }
            else if (((*i).second.second == "") != ((*j).second.second == ""))
            {
                // std::cout << "[3] " << lhs.serial() << "(" << val1 << ") vs " << rhs.serial() << "(" << val2 << ") -> '" << (*i).second.second << "' < '" << (*j).second.second << "'" << std::endl;
                return ((*i).second.second == ""? 1 : -1);
            }

            ++i;
            ++j;
        }

        if (lhs.num_sites() != rhs.num_sites())
        {
            return (lhs.num_sites() < rhs.num_sites()? 1 : -1);
        }

        ignores_.insert(std::lower_bound(ignores_.begin(), ignores_.end(), pair_key), pair_key);
        i = lhs.begin();
        j = rhs.begin();
        while (i != lhs.end() && j != rhs.end())
        {
            if ((*i).second.second != "" && (*i).second.second != "")
            {
                const std::vector<site_type>& pair1(connections_[(*i).second.second]);
                const std::vector<site_type>& pair2(connections_[(*j).second.second]);
                const site_type& target1(
                    (pair1[0].first == val1 && pair1[0].second == (*i).first)? pair1[1] : pair1[0]);
                const site_type& target2(
                    (pair2[0].first == val2 && pair2[0].second == (*j).first)? pair2[1] : pair2[0]);
                if (target1.second != target2.second)
                {
                    ignores_.pop_back();
                    return (target1.second < target2.second? 1 : -1);
                }
                const int retval(compare(target1.first, target2.first));
                // std::cout << "[0] " << lhs.serial() << "(" << val1 << ") vs " << rhs.serial() << "(" << val2 << ") -> " << retval << std::endl;
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

protected:

    const Species& root_;
    connection_container_type connections_;
    std::vector<std::pair<index_type, index_type> > ignores_;
};

std::string serialize_species(const Species& sp)
{
    unit_species_comparerator comp(sp);
    std::vector<unit_species_comparerator::index_type> units;
    for (unit_species_comparerator::index_type i(0); i < sp.num_units(); ++i)
    {
        units.push_back(i);
    }
    std::sort(units.begin(), units.end(), comp);

    std::vector<UnitSpecies> usps(sp.list_units());
    utils::get_mapper_mf<std::string, std::string>::type cache;
    unsigned int stride(1);
    std::stringstream ss;
    for (std::vector<unit_species_comparerator::index_type>::const_iterator i(units.begin());
        i != units.end(); ++i)
    {
        UnitSpecies& usp(usps.at(*i));
        for (UnitSpecies::container_type::size_type j(0); j < usp.num_sites(); ++j)
        {
            UnitSpecies::container_type::value_type& site(usp.at(j));
            if (site.second.second == "" or site.second.second[0] == '_')
            {
                continue;
            }

            utils::get_mapper_mf<std::string, std::string>::type::const_iterator
                it(cache.find(site.second.second));
            if (it == cache.end())
            {
                ss << stride;
                cache.insert(std::make_pair(site.second.second, ss.str()));
                site.second.second = ss.str();
                ++stride;
                ss.clear();
                ss.str("");
            }
            else
            {
                site.second.second = (*it).second;
            }
        }
    }

    Species newsp;
    for (std::vector<UnitSpecies>::const_iterator i(usps.begin()); i != usps.end(); ++i)
    {
        newsp.add_unit(*i);
    }
    return newsp.serial();
}

} // ecell4
