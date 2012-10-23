
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>

#include <vector>
#include <numeric>
#include <map>

#include "Gillespie.hpp"

// UNITTEST:	g++ Gillespie.cpp -lgsl -lgslcbas -DUNITTEST

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
//	 Reaction Class Support Routines
//============================================================
Reaction::Reaction(void) {
	validp = false;
	k = 0.0;
}

void Reaction::check(void) {
	if (0 < substances.size() && 0 < products.size() && k != 0.0) {
		this->validp = true;
	} else {
		this->validp = false;
	}
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
//	GillespieSolver 	*Definitions
//============================================================
GillespieSolver::GillespieSolver(void)
{
	current_t = 0.0;

	// initialize random number generator
	T = gsl_rng_default;
	this->random_handle = gsl_rng_alloc(T);
	gsl_rng_set(this->random_handle, time(NULL));
}

GillespieSolver::~GillespieSolver(void)
{
	// freeing random number handle.
	gsl_rng_free(this->random_handle);
}

// GillespieSolver
//  	*about reactions
int GillespieSolver::reaction_add(void)
{
	// return index of reactions
	int retval = this->models.size();
	Reaction new_react;
	this->models.push_back(new_react);
	return retval;
}

void GillespieSolver::reaction_add_substance(
		int reaction_num,
		int specie_id,
		int stoichiometry )
{
	Reaction *r = &(this->models[reaction_num]);
	if (r == NULL)
		return;
	r->substances.push_back( Specie_index_number(specie_id, stoichiometry) );
	r->check();
}

void GillespieSolver::reaction_add_product(
		int reaction_num,
		int specie_id,
		double stoichiometry )
{
	Reaction *r = &(this->models[reaction_num]);
	if (r == NULL)
		return;
	r->products.push_back( Specie_index_number(specie_id, stoichiometry) );
	r->check();
}

void GillespieSolver::reaction_set_kinetic_parameter(int reaction_num, double k)
{
	Reaction *r = &(this->models[reaction_num]);
	if (r == NULL) 
		return;
	r->k = k;
	r->check();
}

// GillespieSolver
//  	*Properties.
void GillespieSolver::set_current_time(double new_t) 
{	this->current_t = new_t;	}

double GillespieSolver::get_current_time(void) 
{	return this->current_t;		}

int GillespieSolver::get_current_state(int *array, int len)
{
	if (array == NULL)
		return -1;

	int idx = 0;
	std::vector<int>::iterator it = this->current_state.begin();
	while( it != this->current_state.end() && idx < len ) {
		*(array + idx) = *it;
		idx++;
		it++;
	}
	return idx;
}

void GillespieSolver::set_current_state(int *array, int len) {
	this->current_state.clear();
	for(int i = 0; i < len; i++) {
		this->current_state.push_back(*(array + i));
	}
}

// GillespieSolver::step() function returns dt.
double GillespieSolver::step(void)
{
	if (models.size() == 0 || current_state.size() == 0) {
		// reactions or world status not initialized.
		return 0.0;
	}

	std::vector<double>	a(this->models.size() );

	for(size_t idx(0); idx < this->models.size(); idx++) {
		a[idx] = this->models[idx].k;
		for(Species_Vector::iterator it(models[idx].substances.begin());
				it != models[idx].substances.end();
				it++) {
			a[idx] *= combination(
					this->current_state[it->specie_index],
					it->specie_stoichiometry
				);
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

	this->current_t += dt;
	//	Ru(models[u]) occurs.
	for(Species_Vector::iterator it(models[u].substances.begin());
			it != models[u].substances.end();
			it++) {
		this->current_state[it->specie_index] -= it->specie_stoichiometry;
	}
	for(Species_Vector::iterator it(models[u].products.begin());
			it != models[u].products.end();
			it++) {
		this->current_state[it->specie_index] += it->specie_stoichiometry;
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

