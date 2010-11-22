#ifndef STRUCTURE_TYPE_HPP
#define STRUCTURE_TYPE_HPP

#include <ostream>
#include <string>
#include <boost/format.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "exceptions.hpp"
#include "utils/get_mapper_mf.hpp"

class ParticleModel;

class StructureType
{
    friend class ParticleModel;
private:
    typedef get_mapper_mf<std::string, std::string>::type string_map_type;

public:
    typedef std::string identifier_type;
    typedef string_map_type::const_iterator string_map_iterator;
    typedef boost::iterator_range<string_map_iterator> attributes_range;

public:
    identifier_type const& id() const;

    std::string const& operator[](std::string const& name) const;

    std::string& operator[](std::string const& name);

    attributes_range attributes() const;

    ParticleModel* model() const
    {
        return model_;
    }

    StructureType(): model_(0) {}
 
protected:
    void bind_to_model(ParticleModel* model, identifier_type const& id)
    {
        model_ = model; 
        id_ = id;
    }

private:
    ParticleModel* model_;
    identifier_type id_;
    string_map_type attrs_;
};

template<typename Tchar_, typename Ttraits_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& out, const StructureType& v)
{
    bool first = true;
    out << "StructureType(id=" << v.id() << ", attributes={";

    typename StructureType::attributes_range attributes(v.attributes());
    for (typename boost::range_const_iterator<
        typename StructureType::attributes_range>::type
            i(attributes.begin()), e(attributes.end()); i != e; ++i)
    {
        typename boost::range_value<typename StructureType::attributes_range>::type
                const& pair(*i);
        if (!first)
            out << ", ";
        out << pair.first << ":" << pair.second;
        first = false;
    }
    out << "})";
    return out;
}

#endif /* STRUCTURE_TYPE_HPP */
