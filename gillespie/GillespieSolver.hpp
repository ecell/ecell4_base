#include <map>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>


//============================================================
//	World	
//============================================================
class World {
private:
public:	// Data members
	double current_t;
	std::map<int,int> current_state;	// [id] -> number of substrate	// id => 分子の個数　

public:
	World();
	double get_current_time(void);
	void set_current_time(double new_t);
	int get_current_state(int id);
	void set_current_state(int id, int number);
	void add_specie(int id, int number);
};

//============================================================
//	Model and supporting classes.
//============================================================
typedef std::pair<int,int> id_stoichiometry;
class ReactionRule {
private:
	bool valid_react;
	void valid_check(void);

public:	// should be private member?
	std::vector<id_stoichiometry>	reactants;	// id
	std::vector<id_stoichiometry>	products;	// id
	double k;
public:
	void add_reactant(int id, int stoichiometry);
	void add_product(int id, int stoichiometry);
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
