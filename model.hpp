#ifndef MODEL_HPP
#define MODEL_HPP

#include "serial_id_generator.hpp"
#include "species_type_id.hpp"
#include "network_rules.hpp"

class species_type;

class model
{
private:
    typedef serial_id_generator<species_type_id> species_type_id_generator_type;
    typedef std::map<species_type_id, species_type*> species_type_map_type;

public:
    model() {}

    ~model();

    class network_rules& network_rules()
    {
        return network_rules_;
    }

    class network_rules const& network_rules() const
    {
        return network_rules_;
    }

    species_type* new_species_type();

public:
    species_type_id_generator_type species_type_id_generator_;
    species_type_map_type species_type_map_;
    class network_rules network_rules_;
};


#endif /* MODEL_HPP */
