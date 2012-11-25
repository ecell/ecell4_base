
#include <map>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>


#ifndef INCLUDE_GUARD_GILLESPIE_SOLVER
#	define INCLUDE_GUARD_GILLESPIE_SOLVER

using namespace std;

//============================================================
//	Model and supporting classes.
//============================================================
typedef std::pair<std::string,int> id_stoichiometry;
class ReactionRule {
private:
	bool valid_react;
	void valid_check(void);

public:	// should be private member?
	std::vector<id_stoichiometry>	reactants;	// id
	std::vector<id_stoichiometry>	products;	// id
	double k;
public:
	void add_reactant(string &sp, int stoichiometry);
	void add_product(string &sp, int stoichiometry);
	void set_kinetic_parameter(double new_k);
	bool is_valid(void);
};

class Model {
public:		// XXX should be private member ? 
	std::vector<ReactionRule> reactions;
};

//============================================================
//	Gillespie Solver 	*Prototype Declaration.
//============================================================
class GillespieSolver {
private:
	// for random number 
	gsl_rng *random_handle;
	const gsl_rng_type *T;
	Model &m;
	World &w;

public:
	GillespieSolver(World &w, Model &m);
	~GillespieSolver();

	// Functions about reactions.
	double step(void);
	double run(double duration);
};

#endif	//INCLUDE_GUARD_GILLESPIE_SOLVER
