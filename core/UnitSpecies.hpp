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

#include "get_mapper_mf.hpp"


namespace ecell4
{

class UnitSpecies
{
public:

    typedef std::string serial_type;
    typedef std::pair<std::string, std::string> site_type;
    typedef utils::get_mapper_mf<std::string, site_type>::type container_type;

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
        container_type::const_iterator itr(sites_.find(name));
        if (itr != sites_.end())
        {
            return false;
        }
        sites_.insert(std::make_pair(name, std::make_pair(state, bond)));
        return true;
    }

    bool has_site(const std::string& name) const
    {
        return sites_.find(name) != sites_.end();
    }

    const site_type& get_site(const std::string& name) const
    {
        return (*sites_.find(name)).second;
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

