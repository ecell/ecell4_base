#ifndef MODEL_HPP
#define MODEL_HPP

#include "SerialIDGenerator.hpp"
#include "SpeciesTypeID.hpp"
#include "NetworkRules.hpp"

class SpeciesType;

class Model
{
private:
    typedef SerialIDGenerator<SpeciesTypeID> species_type_id_generator_type;
    typedef std::map<SpeciesTypeID, SpeciesType*> species_type_map_type;

public:
    Model();

    ~Model();

    class NetworkRules& network_rules()
    {
        return *network_rules_;
    }

    NetworkRules const& network_rules() const
    {
        return *network_rules_;
    }

    SpeciesType* new_species_type();

public:
    species_type_id_generator_type species_type_id_generator_;
    species_type_map_type species_type_map_;
    NetworkRules *network_rules_;
};


#endif /* MODEL_HPP */
