#ifndef SPECIES_TYPE_HPP
#define SPECIES_TYPE_HPP

#include <ostream>
#include <string>
#include <boost/format.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "SpeciesTypeID.hpp"

#include "exceptions.hpp"
#include "get_mapper_mf.hpp"

class Model;

class SpeciesType
{
    friend class Model;
private:
    typedef get_mapper_mf<std::string, std::string>::type string_map_type;

public:
    typedef string_map_type::const_iterator string_map_iterator;
    typedef boost::iterator_range<string_map_iterator> attributes_range;

public:
    SpeciesTypeID const& id() const
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

    attributes_range attributes() const
    {
        return attributes_range(attrs_.begin(), attrs_.end());
    }

    Model* model() const
    {
        return model_;
    }

    SpeciesType(SpeciesTypeID const& id): id_(id) {}
 
protected:
    class Model*& model()
    {
        return model_;
    }

private:
    class Model* model_;
    SpeciesTypeID id_;
    std::string name_;
    string_map_type attrs_;
};

template<typename Tchar_, typename Ttraits_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& out, const SpeciesType& v)
{
    bool first = true;
    out << "SpeciesType(id=" << v.id() << ", attributes={";

    typename SpeciesType::attributes_range attributes(v.attributes());
    for (typename boost::range_const_iterator<
        typename SpeciesType::attributes_range>::type
            i(attributes.begin()), e(attributes.end()); i != e; ++i)
    {
        typename boost::range_value<typename SpeciesType::attributes_range>::type
                const& pair(*i);
        if (!first)
            out << ", ";
        out << pair.first << ":" << pair.second;
        first = false;
    }
    out << "})";
    return out;
}

#endif /* SPECIES_TYPE_HPP */
