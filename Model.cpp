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
    delete network_rules_;
}

void Model::add_species_type(boost::shared_ptr<species_type_type> const& species)
{
    species->bind_to_model(this, species_type_id_generator_());
    species_type_map_.insert(std::make_pair(species->id(), species));
}

void Model::add_species_type(boost::shared_ptr<species_type_type> const& species);

boost::shared_ptr<Model::species_type_type> Model::get_species_type_by_id(SpeciesTypeID const& id) const
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

std::string const& Model::operator[](std::string const& name) const
{
    string_map_type::const_iterator i(attrs_.find(name));
    if (i == attrs_.end())
        throw not_found((boost::format("key %s not found") % name).str());
    return (*i).second;
}

std::string& Model::operator[](std::string const& name)
{
    return attrs_[name];
}
