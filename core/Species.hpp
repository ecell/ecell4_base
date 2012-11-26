#ifndef __SPECIES_HPP
#define __SPECIES_HPP

#include <string>
#include <vector>
#include <map>

#include "types.hpp"
#include "exceptions.hpp"


namespace ecell4
{

class Species
{
public:

    typedef std::map<std::string, std::string> attributes_container_type;

    Species(std::string const& name = "")
        : name_(name)
    {
        ;
    }

    std::string name() const
    {
        return name_;
    }

    typename attributes_container_type::mapped_type get_attribute(
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
        typename attributes_container_type::mapped_type value)
    {
        attributes_[name_attr] = value;
    }

    bool operator==(Species const& rhs) const;
    bool operator<(Species const& rhs) const;
    bool operator>(Species const& rhs) const;

protected:

    std::string name_;
    attributes_container_type attributes_;
};

typedef std::vector<Species> SpeciesVector;

} // ecell4

#endif /* __SPECIES_HPP */
