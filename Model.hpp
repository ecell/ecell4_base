#ifndef MODEL_HPP
#define MODEL_HPP

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "SerialIDGenerator.hpp"
#include "SpeciesTypeID.hpp"
#include "SpeciesType.hpp"
#include "NetworkRules.hpp"
#include "utils/get_mapper_mf.hpp"
#include "utils/pair.hpp"

#include "SpeciesType.hpp"

class Model: private boost::noncopyable
{
public:
    typedef SpeciesTypeID species_id_type;
    typedef SpeciesType species_type_type;

private:
    typedef SerialIDGenerator<species_id_type> species_type_id_generator_type;
    typedef std::map<species_id_type, boost::shared_ptr<species_type_type> > species_type_map_type;
    typedef select_second<species_type_map_type::value_type> second_selector_type;
    typedef get_mapper_mf<std::string, std::string>::type string_map_type;

public:
    typedef boost::transform_iterator<second_selector_type,
            species_type_map_type::const_iterator> species_type_iterator;
    typedef boost::iterator_range<species_type_iterator> species_type_range;
    typedef NetworkRules network_rules_type;
    typedef string_map_type::const_iterator string_map_iterator;
    typedef boost::iterator_range<string_map_iterator> attributes_range;

public:
    Model();

    ~Model();

    NetworkRules& network_rules() const
    {
        return *network_rules_;
    }

    void add_species_type(boost::shared_ptr<species_type_type> const& species);

    boost::shared_ptr<species_type_type> get_species_type_by_id(species_id_type const& id) const;

    species_type_range get_species_types() const;

    std::string const& operator[](std::string const& name) const;

    std::string& operator[](std::string const& name);

    attributes_range attributes() const
    {
        return attributes_range(attrs_.begin(), attrs_.end());
    }


public:
    species_type_id_generator_type species_type_id_generator_;
    species_type_map_type species_type_map_;
    NetworkRules *network_rules_;
    string_map_type attrs_;
};


#endif /* MODEL_HPP */
