
#define BOOST_TEST_MODULE "GillespieWorld_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <core/RandomNumberGenerator.hpp>
#include <core/Model.hpp>
#include <core/NetworkModel.hpp>

#include "../GillespieWorld.cpp"
#include "../GillespieSimulator.hpp"

using namespace ecell4;
using namespace ecell4::gillespie;

BOOST_AUTO_TEST_CASE(GillespieWorld_test)
{
	GSLRandomNumberGenerator rng;
	rng.seed(time(NULL));
	Species sp1("A");
	Species sp2("B");

	
	Real vol(1.0);
	boost::shared_ptr<GillespieWorld> world(new GillespieWorld(vol));

	world->add_species(sp1);
	world->add_species(sp2);

	world->add_molecules(sp1, 10);
	world->add_molecules(sp2, 20);
	world->set_t(0.5);

	BOOST_CHECK(world->t() == 0.5);
	BOOST_CHECK(world->num_species() == 2);
	BOOST_CHECK(world->num_molecules(sp1) == 10);
	BOOST_CHECK(world->num_molecules(sp2) == 20);

}
