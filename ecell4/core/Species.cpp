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
    throw NotSupported("Function 'Species::count' was deprecated. Rather use 'count_species_matches'");
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

Species::attribute_type Species::get_attribute(const std::string& key) const
{
    return attributes_.get(key);
}

std::vector<std::pair<std::string, Species::attribute_type> > Species::list_attributes() const
{
    return attributes_.values();
}

void Species::set_attributes(const Species& sp)
{
    attributes_ = sp.attributes();
}

void Species::overwrite_attributes(const Species& sp)
{
    attributes_.overwrite(sp.attributes());
}

void Species::remove_attribute(const std::string& key)
{
    attributes_.remove(key);
}

bool Species::has_attribute(const std::string& key) const
{
    return attributes_.has_key(key);
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

} // ecell4
