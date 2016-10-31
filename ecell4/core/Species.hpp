#ifndef __ECELL4_SPECIES_HPP
#define __ECELL4_SPECIES_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

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

protected:

    typedef utils::get_mapper_mf<std::string, std::string>::type
    attributes_container_type;

public:

    Species()
        : units_()
    {
        ; // do nothing
    }

    // Species(const Species& sp)
    //     : units_()
    // {
    //     deserialize(sp.serial());
    // }

    explicit Species(const serial_type& name)
        : units_()
    {
        deserialize(name);
    }

    // Species(
    //     const serial_type& name, const std::string& D)
    //     : units_()
    // {
    //     deserialize(name);
    //     set_attribute("D", D);
    // }

    Species(
        const serial_type& name, const std::string& radius, const std::string& D,
        const std::string location = "")
        : units_()
    {
        deserialize(name);
        set_attribute("radius", radius);
        set_attribute("D", D);
        set_attribute("location", location);
    }

    void deserialize(const serial_type& serial);
    serial_type serial() const;

    Integer num_units() const
    {
        return units_.size();
    }

    void add_unit(const UnitSpecies& usp);

    inline container_type::const_iterator begin() const
    {
        return units_.begin();
    }

    inline container_type::const_iterator end() const
    {
        return units_.end();
    }

    const std::vector<UnitSpecies>& units() const
    {
        return units_;
    }

    const UnitSpecies& at(const container_type::size_type& idx) const
    {
        return units_.at(idx);
    }

    // Integer get_unit(const UnitSpecies& usp)
    // {
    //     container_type::iterator itr;
    //     for (itr = units_.begin(); itr != units_.end(); ++itr)
    //     {
    //         if (usp == *itr)
    //         {
    //             return itr - units_.begin();
    //         }
    //     }
    //     throw NotFound("UnitSpecies not found");
    // }

    // const std::vector<UnitSpecies> list_sites()
    // {
    //     std::vector<UnitSpecies> usps;
    //     if (units_.size() == 0)
    //     {
    //         return usps;
    //     }
    //     container_type::const_iterator it(units_.begin());
    //     ++it;
    //     for (; it != units_.end(); ++it)
    //     {
    //           if ((*it).sites_.size() != 0)
    //           {
    //             usps.push_back((*it).serial());
    //           }
    //     }
    //     return usps;
    // }

    const attributes_container_type& attributes() const
    {
        return attributes_;
    }

    std::vector<std::pair<std::string, std::string> > list_attributes();
    std::string get_attribute(const std::string& name_attr) const;
    void set_attribute(const std::string& name_attr, const std::string& value);
    void set_attributes(const Species& sp);
    void overwrite_attributes(const Species& sp);
    void remove_attribute(const std::string& name_attr);
    bool has_attribute(const std::string& name_attr) const;

    bool operator==(const Species& rhs) const;
    bool operator!=(const Species& rhs) const;
    bool operator<(const Species& rhs) const;
    bool operator>(const Species& rhs) const;

    Integer count(const Species& sp) const;

    /** for epdp
     */
    serial_type name() const
    {
        return serial();
    }

protected:

    std::vector<UnitSpecies> units_;
    attributes_container_type attributes_;
};

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

#endif /* __ECELL4_SPECIES_HPP */
