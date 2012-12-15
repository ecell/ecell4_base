
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
	for(std::vector<double>::iterator it = v.begin(); it != v.end(); it++) 
	{
		if (first == true) 
		{	
			first = false;
			std::cout << "{ ";
		} 
		else 
		{
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
	for(std::vector<int>::iterator it = v.begin(); it != v.end(); it++) 
	{
		if (first == true) 
		{	
			first = false;
			std::cout << "{ ";
		} 
		else 
		{
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
int factorial(int n) 
{
	int product(1);
	for(int i(1); i <= n; i++) 
	{
		product *= i;
	}
	return product;
}

// calcurate nPk
int permutation(int n, int k) 
{
	int ans(1);
	for(int i(0); i < k; i++) 
	{
		ans *= n;
		n--;
	}
	return ans;
}

int combination(int n, int k) 
{
	int kk = k < (n - k) ? k : n-k;
	return permutation(n, kk) / factorial(kk);
}


//============================================================
//	Gillespie Solver 	*New Implementation
//============================================================
namespace ecell4 
{
namespace gillespie 
{

void GillespieSolver::step(void) 
{
	const NetworkModel::reaction_rule_container_type &possible_reaction_rules = this->model_.reaction_rules();

	if (possible_reaction_rules.size() == 0) 
	{
		printf("Error %s::%d\n", __FILE__, __LINE__);
		return;	// Throw an exception?
	}
	std::vector<double> a( possible_reaction_rules.size() ) ;

	for(unsigned int idx(0); idx < possible_reaction_rules.size(); idx++) 
	{
		a[idx] = possible_reaction_rules[idx].k();
		const ReactionRule::reactant_container_type &reactants = possible_reaction_rules[idx].reactants();
		for( ReactionRule::reactant_container_type::iterator it = reactants.begin();
				it != reactants.end();
				it++ )
		{
			// assumption: stoichiometry = 1.
			a[idx] *= this->world_.num_molecules(*it);
		}
	}
	double a_total( std::accumulate(a.begin(), a.end(), double(0.0)) );
	if(a_total == 0.0) 
	{
		return;		// XXX Throw an exception?
	}

	double rnd_num1 = this->rng_.uniform(0, 1);
	double dt = gsl_sf_log(1.0 / rnd_num1) / double(a_total);
	double rnd_num2 = this->rng_.uniform(0, 1) * a_total;

	int u(-1);
	double acc(0.0);
	int len( a.size() );
	do 
	{
		u++;
		acc += a[u];
	} while(acc < rnd_num2 && u < len - 1);

	// Reaction[u] occures.
	for(
		ReactionRule::reactant_container_type::iterator it(possible_reaction_rules[u].reactants().begin());
		it != possible_reaction_rules[u].reactants().end();
		it++ )
	{
		int one(1);
		this->world_.remove_molecules(*it, one);
	}
	for(
		ReactionRule::product_container_type::iterator it(possible_reaction_rules[u].products().begin());
		it != possible_reaction_rules[u].products().end();
		it++ ) 
	{
		int one(1);
		this->world_.add_molecules(*it, one);
	}
	this->world_.set_t( this->world_.t() + dt );
}

void GillespieSolver::run(double duration) 
{
	;
}


}	// gillespie

}	// ecell4
