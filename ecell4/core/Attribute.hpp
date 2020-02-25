#ifndef ECELL4_ATTRIBUTE_HPP
#define ECELL4_ATTRIBUTE_HPP

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/variant.hpp>
#include <boost/container/flat_map.hpp>
// #include <boost/container/small_vector.hpp>

#include <ecell4/core/config.h>

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

    std::vector<value_type> values() const;
    mapped_type get(const key_type& name_attr) const;

    template <typename T_>
    T_ get_as(const key_type& name_attr) const
    {
        mapped_type val = get(name_attr);
        if (T_* x = boost::get<T_>(&val))
        {
            return (*x);
        }
        throw NotSupported("An attribute has incorrect type.");
    }

    template <typename T_>
    void set(const key_type& name_attr, T_ value)
    {
        data_[name_attr] = value;
    }

    void remove(const key_type& name_attr);
    bool has_key(const key_type& name_attr) const;
    void overwrite(const Attribute& attr);
    void clear();

protected:

    container_type data_;
};

template <> Real Attribute::get_as<Real>(const key_type& name_attr) const;
template <> Integer Attribute::get_as<Integer>(const key_type& name_attr) const;
template <> void Attribute::set<const char*>(const key_type& name_attr, const char* value);
template <> void Attribute::set<Real>(const key_type& name_attr, const Real value);
template <> void Attribute::set<Integer>(const key_type& name_attr, const Integer value);

} // ecell4


#endif /* ECELL4_ATTRIBUTE_HPP */
