#include "Species.hpp"
#include "Context.hpp"

#include <algorithm>


namespace ecell4
{

bool Species::operator==(const Species& rhs) const
{
    return (serial() == rhs.serial());
}

bool Species::operator!=(const Species& rhs) const
{
    return (serial() != rhs.serial());
}

bool Species::operator<(const Species& rhs) const
{
    return (serial() < rhs.serial());
}

bool Species::operator>(const Species& rhs) const
{
    return (serial() > rhs.serial());
}

Integer Species::count(const Species& sp) const
{
    return count_spmatches(*this, sp);
}

void Species::add_unit(const UnitSpecies& usp)
{
    if (usp.name() == "")
    {
        throw NotSupported("UnitSpecies must have a name.");
    }
    else if (serial_ != "")
    {
        serial_ += "." + usp.serial();
    }
    else
    {
        serial_ = usp.serial();
    }
}

std::vector<std::pair<std::string, Species::attribute_type> > Species::list_attributes() const
{
    std::vector<std::pair<std::string, attribute_type> > retval;
    for (attributes_container_type::const_iterator
        i(attributes_.begin()); i != attributes_.end(); ++i)
    {
        retval.push_back(*i);
    }
    return retval;
}

Species::attribute_type Species::get_attribute(const std::string& name_attr) const
{
    attributes_container_type::const_iterator
        i(attributes_.find(name_attr));
    if (i == attributes_.end())
    {
        std::ostringstream message;
        message << "attribute [" << name_attr << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }

    return (*i).second;
}

// std::string Species::get_attribute(const std::string& name_attr) const
// {
//     return boost::get<std::string>(get_attribute_as_variant(name_attr));
// }

void Species::set_attributes(const Species& sp)
{
    attributes_ = sp.attributes();
}

void Species::overwrite_attributes(const Species& sp)
{
    const attributes_container_type& attrs(sp.attributes());
    for (attributes_container_type::const_iterator i(attrs.begin());
        i != attrs.end(); ++i)
    {
        this->set_attribute((*i).first, (*i).second);
    }
}

void Species::remove_attribute(const std::string& name_attr)
{
    attributes_container_type::iterator
        i(attributes_.find(name_attr));
    if (i == attributes_.end())
    {
        std::ostringstream message;
        message << "attribute [" << name_attr << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }

    attributes_.erase(i);
}

bool Species::has_attribute(const std::string& name_attr) const
{
    return (attributes_.find(name_attr) != attributes_.end());
}

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
                if ((*i).second.second == "" || is_wildcard((*i).second.second))
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
            if ((*i).second.second == "" || is_wildcard((*i).second.second))
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

Species format_species(const Species& sp)
{
    unit_species_comparerator comp(sp);
    const std::vector<UnitSpecies>::size_type num_units = comp.units().size();

    std::vector<unit_species_comparerator::index_type> units;
    for (unit_species_comparerator::index_type i(0); i < num_units; ++i)
    {
        units.push_back(i);
    }

    std::sort(units.begin(), units.end(), comp);

    std::vector<unit_species_comparerator::index_type>
        next(num_units, num_units);
    unsigned int stride(0);
    for (unit_species_comparerator::index_type i(0); i < num_units; ++i)
    {
        const unit_species_comparerator::index_type idx(units[i]);
        comp.reorder_units(next, idx, stride);
    }
    for (unsigned int i(0); i < num_units; ++i)
    {
        units[next[i]] = i;
    }

    Species newsp;
    utils::get_mapper_mf<std::string, std::string>::type cache;
    stride = 1;
    std::stringstream ss;
    for (std::vector<unit_species_comparerator::index_type>::const_iterator
        i(units.begin()); i != units.end(); ++i)
    {
        UnitSpecies usp(sp.units().at(*i));
        for (UnitSpecies::container_type::size_type j(0);
            j < static_cast<UnitSpecies::container_type::size_type>(usp.num_sites()); ++j)
        {
            UnitSpecies::container_type::value_type& site(usp.at(j));
            if (site.second.second == "" || is_wildcard(site.second.second))
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
        newsp.add_unit(usp);
    }
    return newsp;
}

} // ecell4
