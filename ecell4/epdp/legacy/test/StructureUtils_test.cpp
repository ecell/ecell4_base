#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "StructureUtils"

#include <boost/test/included/unit_test.hpp>
#include "Box.hpp"
#include "Cylinder.hpp"
#include "World.hpp"
#include "EGFRDSimulator.hpp"
#include "StructureUtils.hpp"

BOOST_AUTO_TEST_CASE(test_random_position)
{
    typedef World<CyclicWorldTraits<Real, Real> > world_type;
    typedef EGFRDSimulatorTraitsBase<world_type> simulator_traits_type;
    typedef EGFRDSimulator<simulator_traits_type> simulator_type;
    typedef world_type::position_type position_type;
    typedef StructureUtils<simulator_type> structure_utils_type;

    world_type::traits_type::rng_type rng;

    {
        boost::scoped_ptr<simulator_type::cylindrical_surface_type>
            cyl_surface(
                structure_utils_type::create_cylindrical_surface(
                    "test",
                    create_vector<position_type>(1., 1., 1.),
                    .5,
                    create_vector<position_type>(0., 0., -1.),
                    1.));

        for (int i = 10000; --i >= 0;) {
            position_type p(structure_utils_type::random_position(*cyl_surface, rng));
            BOOST_CHECK(p[0] == 1.);
            BOOST_CHECK(p[1] == 1.);
            BOOST_CHECK(p[2] <  2.);
            BOOST_CHECK(p[2] >= -2.);
        }
    }

    {
        boost::scoped_ptr<simulator_type::cuboidal_region_type>
            cube_surface(
                structure_utils_type::create_cuboidal_region(
                    "test",
                    create_vector<position_type>(0., 0., 0.),
                    array_gen(1., 1., 1.)));

        for (int i = 10000; --i >= 0;) {
            position_type p(structure_utils_type::random_position(*cube_surface, rng));
            BOOST_CHECK(p[0] >= 0.);
            BOOST_CHECK(p[1] >= 0.);
            BOOST_CHECK(p[2] >= 0.);
            BOOST_CHECK(p[0] < 1.);
            BOOST_CHECK(p[1] < 1.);
            BOOST_CHECK(p[2] < 1.);
        }
    }
}
