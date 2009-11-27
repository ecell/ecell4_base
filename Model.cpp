#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/lexical_cast.hpp>

#include "utils/fun_wrappers.hpp"

#include "SpeciesType.hpp"
#include "BasicNetworkRulesImpl.hpp"
#include "Model.hpp"

Model::Model(): network_rules_(new BasicNetworkRulesImpl())
{
}

Model::~Model()
{
    std::for_each(species_type_map_.begin(), species_type_map_.end(),
            compose_unary(
                delete_ptr<boost::remove_pointer<species_type_map_type::mapped_type>::type>(),
                select_second<species_type_map_type::value_type>()));
    delete network_rules_;
}

SpeciesType* Model::new_species_type()
{
    SpeciesType* retval = new SpeciesType(species_type_id_generator_());
    species_type_map_.insert(std::make_pair(retval->id(), retval));
    return retval;
}

SpeciesType* Model::get_species_type_by_id(SpeciesTypeID const& id) const
{
    species_type_map_type::const_iterator i(species_type_map_.find(id));
    if (species_type_map_.end() == i)
    {
        throw not_found(boost::lexical_cast<std::string>(id));
    }

    return (*i).second;
}

Model::species_type_range Model::get_species_types() const
{
    return species_type_range(
        species_type_iterator(species_type_map_.begin(), second_selector_type()),
        species_type_iterator(species_type_map_.end(), second_selector_type()));
}
