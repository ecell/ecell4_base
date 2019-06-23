#include "Species.hpp"

#include <algorithm>


namespace ecell4
{

Species::Species()
    : serial_(""), attributes_()
{
    ; // do nothing
}

Species::Species(const serial_type& name)
    : serial_(name), attributes_()
{
    ;
}

Species::Species(const Species& another)
    : serial_(another.serial()), attributes_(another.attributes_)
{
    ;
}

Species& Species::operator=(const Species& another)
{
    serial_ = another.serial_;
    attributes_ = another.attributes_;
    return *this;
}

Species::Species(
    const serial_type& name, const Real& radius, const Real& D,
    const std::string location, const Integer& dimension)
    : serial_(name), attributes_()
{
    set_attribute("radius", radius);
    set_attribute("D", D);
    set_attribute("location", location);
    set_attribute("dimension", dimension);
}

Species::Species(
    const serial_type& name, const Quantity<Real>& radius, const Quantity<Real>& D,
    const std::string location, const Integer& dimension)
    : serial_(name), attributes_()
{
    set_attribute("radius", radius);
    set_attribute("D", D);
    set_attribute("location", location);
    set_attribute("dimension", dimension);
}

const Species::serial_type Species::serial() const
{
    return serial_;
}

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

const std::vector<UnitSpecies> Species::units() const
{
    std::vector<std::string> unit_serials;
    boost::split(unit_serials, serial_, boost::is_any_of("."));

    std::vector<UnitSpecies> units_;
    for (std::vector<std::string>::const_iterator i(unit_serials.begin());
        i != unit_serials.end(); ++i)
    {
        UnitSpecies usp;
        usp.deserialize(*i);
        units_.insert(std::lower_bound(units_.begin(), units_.end(), usp), usp);
    }
    return units_;
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

const Species::attributes_container_type& Species::attributes() const
{
    return attributes_;
}

Species& Species::D(const std::string& value)
{
    set_attribute("D", value);
    return (*this);
}

Species* Species::D_ptr(const std::string& value)
{
    return &(this->D(value));
}

Species& Species::radius(const std::string& value)
{
    set_attribute("radius", value);
    return (*this);
}

Species* Species::radius_ptr(const std::string& value)
{
    return &(this->radius(value));
}

Species& Species::location(const std::string& value)
{
    set_attribute("location", value);
    return (*this);
}

Species* Species::location_ptr(const std::string& value)
{
    return &(this->location(value));
}

Species& Species::dimension(const std::string& value)
{
    set_attribute("dimension", value);
    return (*this);
}

Species* Species::dimension_ptr(const std::string& value)
{
    return &(this->location(value));
}

/** for epdp
 */
Species::serial_type Species::name() const
{
    return serial();
}

template <>
Real Species::get_attribute_as<Real>(const std::string& name_attr) const
{
    attribute_type val = get_attribute(name_attr);
    if (Quantity<Real>* x = boost::get<Quantity<Real> >(&val))
    {
        return (*x).magnitude;
    }
    else if (Quantity<Integer>* x = boost::get<Quantity<Integer> >(&val))
    {
        return static_cast<Real>((*x).magnitude);
    }
    else if (std::string* x = boost::get<std::string>(&val))
    {
        return std::atof((*x).c_str());
    }
    throw NotSupported("An attribute has incorrect type. Real is expected");
}

template <>
Integer Species::get_attribute_as<Integer>(const std::string& name_attr) const
{
    attribute_type val = get_attribute(name_attr);
    if (Quantity<Integer>* x = boost::get<Quantity<Integer> >(&val))
    {
        return (*x).magnitude;
    }
    else if (std::string* x = boost::get<std::string>(&val))
    {
        return std::atoi((*x).c_str());
    }
    throw NotSupported("An attribute has incorrect type. Integer is expected");
}

template <>
void Species::set_attribute<const char*>(const std::string& name_attr, const char* value)
{
    set_attribute(name_attr, std::string(value));
}

template <>
void Species::set_attribute<Real>(const std::string& name_attr, const Real value)
{
    set_attribute(name_attr, Quantity<Real>(value));
}

template <>
void Species::set_attribute<Integer>(const std::string& name_attr, const Integer value)
{
    set_attribute(name_attr, Quantity<Integer>(value));
}

} // ecell4
