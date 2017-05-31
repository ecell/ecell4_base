#ifndef ECELL4_SPECIES_HPP
#define ECELL4_SPECIES_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/variant.hpp>

#include <ecell4/core/config.h>

#include "hash.hpp"
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

    typedef boost::variant<std::string, Real, Integer, bool> attribute_type;

protected:

    typedef utils::get_mapper_mf<std::string, attribute_type>::type
        attributes_container_type;

public:

    Species()
        : serial_(""), attributes_()
    {
        ; // do nothing
    }

    explicit Species(const serial_type& name)
        : serial_(name), attributes_()
    {
        ;
    }

    Species(const Species& another)
        : serial_(another.serial()), attributes_()
    {
        const std::vector<std::pair<std::string, attribute_type> > attrs = another.list_attributes();
        for (std::vector<std::pair<std::string, attribute_type> >::const_iterator
            i(attrs.begin()); i != attrs.end(); i++)
        {
            set_attribute((*i).first, (*i).second);
        }
    }

    Species(
        const serial_type& name, const Real& radius, const Real& D,
        const std::string location = "")
        : serial_(name), attributes_()
    {
        set_attribute("radius", radius);
        set_attribute("D", D);
        set_attribute("location", location);
    }

    /*
     * The following constructor will be deprecated. Use the above one.
     */
    Species(
        const serial_type& name, const std::string& radius, const std::string& D,
        const std::string location = "")
        : serial_(name), attributes_()
    {
        set_attribute("radius", radius);
        set_attribute("D", D);
        set_attribute("location", location);
    }

    const serial_type serial() const
    {
        return serial_;
    }

    void add_unit(const UnitSpecies& usp);

    const std::vector<UnitSpecies> units() const
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

    const attributes_container_type& attributes() const
    {
        return attributes_;
    }

    std::vector<std::pair<std::string, attribute_type> > list_attributes() const;

    attribute_type get_attribute(const std::string& name_attr) const;

    template <typename T_>
    T_ get_attribute_as(const std::string& name_attr) const
    {
        attribute_type val = get_attribute(name_attr);
        if (T_* x = boost::get<T_>(&val))
        {
            return (*x);
        }
        throw NotSupported("An attribute has incorrect type.");
    }


    template <typename T_>
    void set_attribute(const std::string& name_attr, T_ value)
    // void set_attribute(const std::string& name_attr, const T_& value)
    {
        attributes_[name_attr] = value;
    }

    void set_attributes(const Species& sp);
    void remove_attribute(const std::string& name_attr);
    bool has_attribute(const std::string& name_attr) const;
    void overwrite_attributes(const Species& sp);

    bool operator==(const Species& rhs) const;
    bool operator!=(const Species& rhs) const;
    bool operator<(const Species& rhs) const;
    bool operator>(const Species& rhs) const;

    Integer count(const Species& sp) const;

    /** Method chaining
     */

    Species& D(const std::string& value)
    {
        set_attribute("D", value);
        return (*this);
    }

    inline Species* D_ptr(const std::string& value)
    {
        return &(this->D(value));
    }

    Species& radius(const std::string& value)
    {
        set_attribute("radius", value);
        return (*this);
    }

    inline Species* radius_ptr(const std::string& value)
    {
        return &(this->radius(value));
    }

    Species& location(const std::string& value)
    {
        set_attribute("location", value);
        return (*this);
    }

    inline Species* location_ptr(const std::string& value)
    {
        return &(this->location(value));
    }

    /** for epdp
     */
    serial_type name() const
    {
        return serial();
    }

protected:

    serial_type serial_;
    attributes_container_type attributes_;
};

template <>
inline Real Species::get_attribute_as<Real>(const std::string& name_attr) const
{
    attribute_type val = get_attribute(name_attr);
    if (Real* x = boost::get<Real>(&val))
    {
        return (*x);
    }
    else if (Integer* x = boost::get<Integer>(&val))
    {
        return static_cast<Real>(*x);
    }
    else if (std::string* x = boost::get<std::string>(&val))
    {
        return std::atof((*x).c_str());
    }
    throw NotSupported("An attribute has incorrect type. Real is expected");
}

template <>
inline void Species::set_attribute(const std::string& name_attr, const char* value)
{
    attributes_[name_attr] = std::string(value);
}

Species format_species(const Species& sp);

inline Species::serial_type unique_serial(const Species& sp)
{
    return format_species(sp).serial();
}

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm,
    const ecell4::Species& sp)
{
    strm << sp.serial();
    return strm;
}

} // ecell4

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<ecell4::Species>
{
    std::size_t operator()(const ecell4::Species& val) const
    {
        return hash<ecell4::Species::serial_type>()(val.serial());
    }
};

ECELL4_DEFINE_HASH_END()

#endif /* ECELL4_SPECIES_HPP */
