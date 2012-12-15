
#include <map>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>



#ifndef INCLUDE_GUARD_GILLESPIE_SOLVER
#	define INCLUDE_GUARD_GILLESPIE_SOLVER

#include <core/CompartmentSpace.hpp>
#include <core/NetworkModel.hpp>
#include <core/Species.hpp>
#include "./GillespieSimulator.hpp"

using namespace std;

namespace ecell4 {

namespace gillespie {

class GillespieSolver {
public:
	GillespieSolver(
		NetworkModel &model, GillespieWorld &world, RandomNumberGenerator &rng)
		: model_(model), world_(world), rng_(rng)
	{;}
	void step(void);
	void run(double duration);

protected:
	NetworkModel &model_;
	GillespieWorld &world_;
	RandomNumberGenerator &rng_;
};

}	// gillespie

}	// ecell4

#endif	//INCLUDE_GUARD_GILLESPIE_SOLVER
