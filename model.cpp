#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <boost/type_traits/remove_pointer.hpp>

#include "utils.hpp"

#include "species_type.hpp"
#include "basic_network_rules_impl.hpp"
#include "model.hpp"

model::model(): network_rules_(new basic_network_rules_impl())
{
}

model::~model()
{
    std::for_each(species_type_map_.begin(), species_type_map_.end(),
            compose_unary(
                delete_ptr<boost::remove_pointer<species_type_map_type::mapped_type>::type>(),
                select_second<species_type_map_type::value_type>()));
    delete network_rules_;
}

species_type* model::new_species_type()
{
    species_type* retval = new species_type(species_type_id_generator_());
    species_type_map_.insert(std::make_pair(retval->id(), retval));
    return retval;
}

