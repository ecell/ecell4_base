#ifndef __ECELL4_SPECIES_HPP
#define __ECELL4_SPECIES_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "config.h"

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "UnitSpecies.hpp"


namespace ecell4
{

class Species
{
public:

    typedef UnitSpecies::serial_type serial_type; //XXX: std::string
    typedef std::vector<UnitSpecies> container_type;

protected:

    typedef utils::get_mapper_mf<std::string, std::string>::type
    attributes_container_type;

public:

    Species()
        : units_()
    {
        ; // do nothing
    }

    Species(const serial_type& name)
        : units_()
    {
        deserialize(name);
    }

    Species(
        const serial_type& name, const std::string& D)
        : units_()
    {
        deserialize(name);
        set_attribute("D", D);
    }

    Species(
        const serial_type& name, const std::string& radius, const std::string& D)
        : units_()
    {
        deserialize(name);
        set_attribute("radius", radius);
        set_attribute("D", D);
    }

    void deserialize(const serial_type& serial)
    {
        std::vector<std::string> unit_serials;
        boost::split(unit_serials, serial, boost::is_any_of("."));

        units_.clear();
        for (std::vector<std::string>::const_iterator i(unit_serials.begin());
            i != unit_serials.end(); ++i)
        {
            UnitSpecies usp;
            usp.deserialize(*i);
            add_unit(usp);
        }
    }

    serial_type serial() const
    {
        if (units_.size() == 0)
        {
            return "";
        }

        container_type::const_iterator it(units_.begin());
        serial_type retval((*it).serial());
        ++it;
        for (; it != units_.end(); ++it)
        {
            retval += ".";
            retval += (*it).serial();
        }
        return retval;
    }

    Integer num_units() const
    {
        return units_.size();
    }

    void add_unit(const UnitSpecies& usp)
    {
        if (usp.name() == "")
        {
            throw NotSupported("UnitSpecies must have a name.");
        }
        units_.insert(std::lower_bound(units_.begin(), units_.end(), usp), usp);
    }

    inline container_type::const_iterator begin() const
    {
        return units_.begin();
    }

    inline container_type::const_iterator end() const
    {
        return units_.end();
    }

    Integer get_unit(const UnitSpecies& usp)
    {
        container_type::iterator itr;
        for (itr = units_.begin(); itr != units_.end(); ++itr)
        {
            if (usp == *itr)
            {
                return itr - units_.begin();
            }
        }
        throw NotFound("UnitSpecies not found");
    }

    const std::vector<UnitSpecies> list_sites()
    {
        std::vector<UnitSpecies> usps;
        if (units_.size() == 0)
        {
            return usps;
        }
        container_type::const_iterator it(units_.begin());
        ++it;
        for (; it != units_.end(); ++it)
        {
//            if ((*it).sites_.size() != 0)
//            {
                usps.push_back((*it).serial());
//            }
        }
        return usps;
    }

    bool match(const Species& target) const
    {
        container_type::const_iterator i(units_.begin()), j(target.begin());
        while (i != units_.end())
        {
            const UnitSpecies& usp(*i);

            j = std::lower_bound(j, target.end(), usp);
            if (j == target.end())
            {
                return false;
            }

            const container_type::const_iterator
                nexti(std::upper_bound(i, units_.end(), usp)),
                nextj(std::upper_bound(j, target.end(), usp));
            if (nextj - j < nexti - i)
            {
                return false;
            }

            i = nexti;
            j = nextj;
        }
        return true;

        // for (container_type::const_iterator i(units_.begin()); i != units_.end(); ++i)
        // {
        //     if (std::count(target.begin(), target.end(), (*i))
        //         < std::count(units_.begin(), units_.end(), (*i))) //XXX:
        //     {
        //         return false;
        //     }
        // }
        // return true;
    }

    const attributes_container_type& attributes() const
    {
        return attributes_;
    }

    std::vector<std::pair<std::string, std::string> > list_attributes()
    {
        std::vector<std::pair<std::string, std::string> > retval;
        for (attributes_container_type::const_iterator
            i(attributes_.begin()); i != attributes_.end(); ++i)
        {
            retval.push_back(*i);
        }
        return retval;
    }

    std::string get_attribute(const std::string& name_attr) const
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

    void set_attribute(const std::string& name_attr, const std::string& value)
    {
        attributes_[name_attr] = value;
    }

    void set_attributes(const Species& sp)
    {
        attributes_ = sp.attributes();
    }

    void remove_attribute(const std::string& name_attr)
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

    bool has_attribute(const std::string& name_attr) const
    {
        return (attributes_.find(name_attr) != attributes_.end());
    }

    bool operator==(const Species& rhs) const;
    bool operator<(const Species& rhs) const;
    bool operator>(const Species& rhs) const;

protected:

    std::vector<UnitSpecies> units_;
    attributes_container_type attributes_;
};

} // ecell4

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std
{

namespace tr1
{
#elif defined(HAVE_STD_HASH)
namespace std
{
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost
{
#endif

template<>
struct hash<ecell4::Species>
{
    std::size_t operator()(const ecell4::Species& val) const
    {
        return hash<ecell4::Species::serial_type>()(val.serial());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} // tr1

} // std
#elif defined(HAVE_STD_HASH)
} // std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // boost
#endif

#endif /* __ECELL4_SPECIES_HPP */
