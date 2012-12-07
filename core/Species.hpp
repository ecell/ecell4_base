#ifndef __SPECIES_HPP
#define __SPECIES_HPP

#include <string>
#include <vector>
#include <map>

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


namespace ecell4
{

class Species
{
public:

    typedef std::string serial_type;
    typedef utils::get_mapper_mf<std::string, std::string>::type
    attributes_container_type;

    Species(std::string const& name = "")
        : name_(name)
    {
        ;
    }

    serial_type serial() const
    {
        return name();
    }

    std::string name() const
    {
        return name_;
    }

    attributes_container_type::mapped_type get_attribute(
        std::string const& name_attr) const
    {
        attributes_container_type::const_iterator
            i(attributes_.find(name_attr));
        if (i == attributes_.end())
        {
            throw NotFound("attribute not found");
        }

        return (*i).second;
    }

    void set_attribute(
        std::string const& name_attr,
        attributes_container_type::mapped_type value)
    {
        attributes_[name_attr] = value;
    }

    void remove_attribute(std::string const& name_attr)
    {
        attributes_container_type::iterator
            i(attributes_.find(name_attr));
        if (i == attributes_.end())
        {
            throw NotFound("attribute not found");
        }

        attributes_.erase(i);
    }

    bool operator==(Species const& rhs) const;
    bool operator<(Species const& rhs) const;
    bool operator>(Species const& rhs) const;

protected:

    std::string name_;
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
    std::size_t operator()(ecell4::Species const& val) const
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

#endif /* __SPECIES_HPP */
