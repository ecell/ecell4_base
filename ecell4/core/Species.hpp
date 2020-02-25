#ifndef ECELL4_SPECIES_HPP
#define ECELL4_SPECIES_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <functional>
#include <boost/algorithm/string.hpp>
#include <boost/variant.hpp>
#include <boost/container/flat_map.hpp>
// #include <boost/container/small_vector.hpp>

#include <ecell4/core/config.h>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "UnitSpecies.hpp"
// #include "Context.hpp"
#include "Quantity.hpp"
#include "Attribute.hpp"


namespace ecell4
{

class Species
{
public:

    typedef UnitSpecies::serial_type serial_type; //XXX: std::string
    typedef std::vector<UnitSpecies> container_type;

public:

    typedef Attribute::mapped_type attribute_type;

public:

    Species();
    explicit Species(const serial_type& name);
    Species(const Species& another);
    Species& operator=(const Species& another);
    Species(const serial_type& name, const Real& radius, const Real& D,
            const std::string location = "", const Integer& dimension = 0);
    Species(const serial_type& name, const Quantity<Real>& radius, const Quantity<Real>& D,
            const std::string location = "", const Integer& dimension = 0);

    const serial_type serial() const;

    void add_unit(const UnitSpecies& usp);
    const std::vector<UnitSpecies> units() const;

    const Attribute& attributes() const;

    std::vector<std::pair<std::string, attribute_type> > list_attributes() const;
    attribute_type get_attribute(const std::string& key) const;

    template <typename T_>
    T_ get_attribute_as(const std::string& key) const
    {
        return attributes_.get_as<T_>(key);
    }

    template <typename T_>
    void set_attribute(const std::string& key, T_ value)
    {
        attributes_.set<T_>(key, value);
    }

    void set_attributes(const Attribute& attributes);
    void set_attributes(const Species& sp);
    void remove_attribute(const std::string& key);
    bool has_attribute(const std::string& key) const;
    void overwrite_attributes(const Species& sp);

    bool operator==(const Species& rhs) const;
    bool operator!=(const Species& rhs) const;
    bool operator<(const Species& rhs) const;
    bool operator>(const Species& rhs) const;

    Integer count(const Species& sp) const;

    /** Method chaining
     */

    Species& D(const std::string& value);
    Species* D_ptr(const std::string& value);
    Species& radius(const std::string& value);
    Species* radius_ptr(const std::string& value);
    Species& location(const std::string& value);
    Species* location_ptr(const std::string& value);
    Species& dimension(const std::string& value);
    Species* dimension_ptr(const std::string& value);

    /** for epdp
     */
    serial_type name() const;

protected:

    serial_type serial_;
    Attribute attributes_;
};

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm,
    const ecell4::Species& sp)
{
    strm << sp.serial();
    return strm;
}

} // ecell4

namespace std {
template<>
struct hash<ecell4::Species>
{
    std::size_t operator()(const ecell4::Species& val) const
    {
        return hash<ecell4::Species::serial_type>()(val.serial());
    }
};
} // std

#endif /* ECELL4_SPECIES_HPP */
