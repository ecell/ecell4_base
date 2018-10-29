#include "Species.hpp"

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
    // return count_spmatches(*this, sp);
    throw NotSupported("Function 'Species::count' is deprecated. Rather use 'count_species_matches'");
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

} // ecell4
