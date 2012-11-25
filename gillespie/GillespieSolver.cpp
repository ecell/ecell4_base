
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>

#include <vector>
#include <numeric>
#include <map>

#include "GillespieWorld.hpp"
#include "GillespieSolver.hpp"


//============================================================
//	Debugging Utility( For Debugger call )
//============================================================
void display_vector_double(std::vector<double> &v)
{
	bool first = true;
	for(std::vector<double>::iterator it = v.begin(); it != v.end(); it++) {
		if (first == true) {	
			first = false;
			std::cout << "{ ";
		} else {
			std::cout << ", ";	
		}
		std::cout << *it;
	}
	std::cout << " }";
	std::cout << std::endl;
}
void display_vector_int(std::vector<int> &v)
{
	bool first = true;
	for(std::vector<int>::iterator it = v.begin(); it != v.end(); it++) {
		if (first == true) {	
			first = false;
			std::cout << "{ ";
		} else {
			std::cout << ", ";	
		}
		std::cout << *it;
	}
	std::cout << " }";
	std::cout << std::endl;
}


//============================================================
//	Math Utility Functions
//============================================================
int factorial(int n) {
	int product(1);
	for(int i(1); i <= n; i++) {
		product *= i;
	}
	return product;
}

// calcurate nPk
int permutation(int n, int k) {
	int ans(1);
	for(int i(0); i < k; i++) {
		ans *= n;
		n--;
	}
	return ans;
}

int combination(int n, int k) {
	int kk = k < (n - k) ? k : n-k;
	return permutation(n, kk) / factorial(kk);
}


//============================================================
//	Model 	*Definitions
//============================================================
void ReactionRule::valid_check(void)
{
	if ( 0 < this->reactants.size() && 0 < this->products.size() && k != 0.0 ) {
		this->valid_react = true;
	} else {
		this->valid_react = false;
	}
}
bool ReactionRule::is_valid(void)
{	return this->valid_react;	}

void ReactionRule::add_reactant(string &sp, int stoichiometry)
{
	this->reactants.push_back(std::pair<string,int>(sp, stoichiometry));
	this->valid_check();
}

void ReactionRule::add_product(string &sp, int stoichiometry)
{
	this->products.push_back(std::pair<string,int>(sp, stoichiometry));
	this->valid_check();
}
void ReactionRule::set_kinetic_parameter(double new_k)
{
	this->k = new_k;
	this->valid_check();
}


//============================================================
//	GillespieSolver 	*Definitions
//============================================================
GillespieSolver::GillespieSolver(World &arg_world, Model &arg_model)
		:m(arg_model), w(arg_world)
{
	T = gsl_rng_default;
	this->random_handle = gsl_rng_alloc(T);
	gsl_rng_set(this->random_handle, time(NULL));
}

GillespieSolver::~GillespieSolver(void)
{
	// freeing random number handle.
	gsl_rng_free(this->random_handle);
}


// GillespieSolver::step() function returns dt.
double GillespieSolver::step(void)
{
	if (this->m.reactions.size() == 0 || this->w.current_state.size() == 0) {
		// reactions or world status not initialized.
		return 0.0;
	}

	std::vector<double>	a( this->m.reactions.size() );

	for(unsigned int idx(0); idx < this->m.reactions.size(); idx++) {
		a[idx] = this->m.reactions[idx].k;	// implement and fix accessor 
		for(
			std::vector<id_stoichiometry>::iterator it_reactant(this->m.reactions[idx].reactants.begin());
			it_reactant != this->m.reactions[idx].reactants.end();
			it_reactant++
		) 
		{
			a[idx] *= combination(
						this->w.current_state[ it_reactant->first ],
						it_reactant->second	);
		}
	}

	double a_total = std::accumulate(a.begin(), a.end(), double(0.0) );

	if (a_total == 0.0) {
		// There are no reactions to heppen.
		return 0.0;
	}

	double rnd_num1 = gsl_rng_uniform(this->random_handle);
	double dt = gsl_sf_log(1.0 / rnd_num1) / double(a_total);
	double rnd_num2 = gsl_rng_uniform(this->random_handle) * a_total;

	int u(-1);
	double acc(0.0);
	int len = a.size();
	do {
		u++;
		acc += a[u];
	} while ( acc < rnd_num2 && u < len - 1);

	this->w.current_t += dt;
	//	Ru(this->m.rections[u]) occurs.
	for(
		std::vector<id_stoichiometry>::iterator it(this->m.reactions[u].reactants.begin());
		it != this->m.reactions[u].reactants.end();
		it++
	   )
	{
		this->w.current_state[it->first] -= it->second;	//second is stoichiomety
	}

	for(
		std::vector<id_stoichiometry>::iterator it(this->m.reactions[u].products.begin());
		it != this->m.reactions[u].products.end();
		it++
	   )
	{
		this->w.current_state[it->first] += it->second;
	}

	return dt;
}

double GillespieSolver::run(double duration) {
	double t_advanced(0.0);
	double step_dt(0.0);
	do {
		step_dt = this->step();
		if (step_dt == 0.00) {
			break;
		}
		t_advanced += step_dt;
	} while (t_advanced < duration);
	return t_advanced;
}

