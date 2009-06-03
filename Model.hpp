#ifndef MODEL_HPP
#define MODEL_HPP

#include <boost/noncopyable.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "SerialIDGenerator.hpp"
#include "SpeciesTypeID.hpp"
#include "NetworkRules.hpp"

class SpeciesType;

class Model: private boost::noncopyable
{
private:
    typedef SerialIDGenerator<SpeciesTypeID> species_type_id_generator_type;
    typedef std::map<SpeciesTypeID, SpeciesType*> species_type_map_type;
    typedef select_second<species_type_map_type::value_type> second_selector_type;

public:
    typedef boost::transform_iterator<second_selector_type,
            species_type_map_type::const_iterator> species_type_iterator;
    typedef boost::iterator_range<species_type_iterator> species_type_range;

public:
    Model();

    ~Model();

    NetworkRules& network_rules() const
    {
        return *network_rules_;
    }

    SpeciesType* new_species_type();

    SpeciesType* get_species_type_by_id(SpeciesTypeID const& id) const;

    species_type_range get_species_types() const;

public:
    species_type_id_generator_type species_type_id_generator_;
    species_type_map_type species_type_map_;
    NetworkRules *network_rules_;
};


#endif /* MODEL_HPP */
