#ifndef __ECELL4_UNIT_SPECIES_HPP
#define __ECELL4_UNIT_SPECIES_HPP

#include <iostream>
#include <string>

#include "config.h"

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "types.hpp"
#include "get_mapper_mf.hpp"

#include <vector>
#include <algorithm>


namespace ecell4
{

class UnitSpecies
{
public:

    typedef std::string serial_type;
    typedef std::pair<std::string, std::string> site_type;
    typedef std::vector<std::pair<std::string, site_type> > container_type;

protected:

    typedef struct
    {
        typedef container_type::value_type value_type;

        bool operator()(const value_type& val1, const value_type& val2)
        {
            return val1.first < val2.first;
        }
    } site_comparerator;

public:

    UnitSpecies(const std::string& name = "")
        : name_(name)
    {
        ;
    }

    std::string name() const
    {
        return name_;
    }

    void deserialize(const serial_type& serial);

    serial_type serial() const;

    void clear();

    bool add_site(const std::string& name,
        const std::string& state, const std::string& bond)
    {
        std::pair<std::string, site_type> val(
            std::make_pair(name, std::make_pair(state, bond)));
        if (std::binary_search(
            sites_.begin(), sites_.end(), val, site_comparerator()))
        {
            return false;
        }
        sites_.insert(
            std::lower_bound(sites_.begin(), sites_.end(), val,
                site_comparerator()),
            val);
        return true;
    }

    Integer num_sites() const
    {
        return sites_.size();
    }

    bool has_site(const std::string& name) const
    {
        return std::binary_search(sites_.begin(), sites_.end(),
            std::make_pair(name, site_type()), site_comparerator());
    }

    const site_type& get_site(const std::string& name) const
    {
        return (*std::lower_bound(sites_.begin(), sites_.end(),
            std::make_pair(name, site_type()), site_comparerator())).second;
    }

    inline container_type::const_iterator begin() const
    {
        return sites_.begin();
    }

    inline container_type::const_iterator end() const
    {
        return sites_.end();
    }

    bool operator==(const UnitSpecies& rhs) const
    {
        return (serial() == rhs.serial());
    }

    bool operator<(const UnitSpecies& rhs) const
    {
        return (serial() < rhs.serial());
    }

    bool operator>(const UnitSpecies& rhs) const
    {
        return (serial() > rhs.serial());
    }

    container_type::value_type& at(const container_type::size_type& idx)
    {
        return sites_.at(idx);
    }

protected:

    std::string name_;
    container_type sites_;
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
struct hash<ecell4::UnitSpecies>
{
    std::size_t operator()(const ecell4::UnitSpecies& val) const
    {
        return hash<ecell4::UnitSpecies::serial_type>()(val.serial());
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

