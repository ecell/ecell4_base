#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "model"

#include <boost/test/included/unit_test.hpp>
#include "model.hpp"
#include "array_helper.hpp"
#include "species_type.hpp"

BOOST_AUTO_TEST_CASE(basic)
{
    model m;

    species_type* s1 = m.new_species_type();
    (*s1)["name"] = "S";
    (*s1)["D"] = "1.5e-12";
    (*s1)["radius"] = "5e-9";

    species_type* s2 = m.new_species_type();
    (*s2)["name"] = "P";
    (*s2)["D"] = "1e-12";
    (*s2)["radius"] = "7e-9";

    m.network_rules().add_reaction_rule(
        new_reaction_rule(array_gen(s1), array_gen(s2), .2));
}
