#ifndef SPECIES_TYPE_HPP
#define SPECIES_TYPE_HPP

#include <ostream>
#include <string>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include <species_type_id.hpp>

#include "exceptions.hpp"
#include "get_mapper_mf.hpp"

class model;

class species_type
{
    friend class model;
private:
    typedef get_mapper_mf<std::string, std::string>::type string_map_type;

public:
    typedef string_map_type::const_iterator string_map_iterator;
    typedef boost::iterator_range<string_map_iterator> attributes_range;

public:
    species_type_id const& id() const
    {
        return id_;
    }

    std::string const& operator[](std::string const& name) const
    {
        string_map_type::const_iterator i(attrs_.find(name));
        if (i == attrs_.end())
            throw not_found((boost::format("key %s not found") % name).str());
    }

    std::string& operator[](std::string const& name)
    {
        return attrs_[name];
    }

    attributes_range attributes()
    {
        return attributes_range(attrs_.begin(), attrs_.end());
    }

    class model* model() const
    {
        return model_;
    }

    species_type(species_type_id const& id): id_(id) {}
 
protected:
    class model*& model()
    {
        return model_;
    }

private:
    class model* model_;
    species_type_id id_;
    std::string name_;
    string_map_type attrs_;
};

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const species_type& v)
{
    strm << "species_type(id=" << v.id() << ", attributes={";
    BOOST_FOREACH(species_type::string_map_iterator::value_type pair,
            v.attributes())
    {
        strm << pair.first << ":" << pair.second << ", ";
    }
    strm << "})";
    return strm;
}

#endif /* SPECIES_TYPE_HPP */
