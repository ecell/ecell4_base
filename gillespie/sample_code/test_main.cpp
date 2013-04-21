#include <string>
#include <sstream>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <cstring>

#include <hdf5.h>
#include <H5Cpp.h>

#include "../GillespieSimulator.hpp"

#include "../../H5Save.hpp"

using namespace ecell4;
using namespace ecell4::gillespie;


using namespace H5;


int main(int argc, char **argv)
{
	GSLRandomNumberGenerator rng;
	rng.seed(time(NULL));

	Species sp1("A");
	Species sp2("B");

	// Expand Reaction Rule.
	ReactionRule rr1;
	rr1.set_k(5.001);
	rr1.add_reactant(sp1);
	rr1.add_product(sp2);

	boost::shared_ptr<NetworkModel> model(new NetworkModel());
	model->add_reaction_rule(rr1);
	Real vol(1.0);

	boost::shared_ptr<GillespieWorld> world(new GillespieWorld(vol));
	// registration
	world->add_species(sp1);
	world->add_species(sp2);

	// set molecule nums.
	world->add_molecules(sp1, 10);
	world->add_molecules(sp2, 10);
	
	model->add_species(sp1);
	model->add_species(sp2);
	
	GillespieSimulator sim(model, world, rng);
	//ecell4_hdf5_manager<GillespieWorld, int>	hdf("GillespieWorld_test.hdf5", model, world, "CompartmentSpace");
	sim.save_hdf5_init("GillespieTest.hdf5");

	for(int i = 0; i < 10; i++) {
		sim.step();

		std::ostringstream ost;
		ost << sim.t();
		printf("t = %f A: %llu, B: %llu \n", sim.t(), world->num_molecules(sp1), world->num_molecules(sp2));
		sim.save_hdf5();
	}

}
