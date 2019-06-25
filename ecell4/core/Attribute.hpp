#ifndef ECELL4_ATTRIBUTE_HPP
#define ECELL4_ATTRIBUTE_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/variant.hpp>
#include <boost/container/flat_map.hpp>
// #include <boost/container/small_vector.hpp>

#include <ecell4/core/config.h>

#include "hash.hpp"
#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "Quantity.hpp"


namespace ecell4
{

class Attribute
{
public:

    typedef std::string key_type;
    typedef boost::variant<std::string, Quantity<Real>, Quantity<Integer>, bool> mapped_type;

protected:

    typedef boost::container::flat_map<key_type, mapped_type>
        container_type;

    // typedef boost::container::small_vector<
    //     std::pair<key_type, mapped_type>, 3
    //         > flat_map_backend_type;
    // typedef boost::container::flat_map<
    //     key_type, mapped_type, std::less<key_type>, flat_map_backend_type
    //         > container_type;

public:

    typedef container_type::value_type value_type;

public:

    Attribute();
    Attribute(const Attribute& another);
    Attribute& operator=(const Attribute& another);
    Attribute(const container_type& attr);

    std::vector<value_type> list_attributes() const;
    mapped_type get_attribute(const key_type& name_attr) const;

    template <typename T_>
    T_ get_attribute_as(const key_type& name_attr) const
    {
        mapped_type val = get_attribute(name_attr);
        if (T_* x = boost::get<T_>(&val))
        {
            return (*x);
        }
        throw NotSupported("An attribute has incorrect type.");
    }

    template <typename T_>
    void set_attribute(const key_type& name_attr, T_ value)
    {
        attributes_[name_attr] = value;
    }

    void set_attributes(const Attribute& attr);
    void remove_attribute(const key_type& name_attr);
    bool has_key(const key_type& name_attr) const;
    void overwrite_attributes(const Attribute& attr);

protected:

    container_type attributes_;
};

template <> Real Attribute::get_attribute_as<Real>(const key_type& name_attr) const;
template <> Integer Attribute::get_attribute_as<Integer>(const key_type& name_attr) const;
template <> void Attribute::set_attribute<const char*>(const key_type& name_attr, const char* value);
template <> void Attribute::set_attribute<Real>(const key_type& name_attr, const Real value);
template <> void Attribute::set_attribute<Integer>(const key_type& name_attr, const Integer value);

} // ecell4


#endif /* ECELL4_ATTRIBUTE_HPP */
