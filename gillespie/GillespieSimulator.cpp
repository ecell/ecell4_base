#include "GillespieSimulator.hpp"
#include <numeric>
#include <vector>
#include <gsl/gsl_sf_log.h>

namespace ecell4 
{

namespace gillespie 
{

void GillespieSimulator::step(void) 
{
	const NetworkModel::reaction_rule_container_type &possible_reaction_rules =
		this->model_->reaction_rules();
	if (possible_reaction_rules.size() == 0)
	{
		return;
	}
	std::vector<double> a(possible_reaction_rules.size() );
	for(unsigned int idx(0); idx < possible_reaction_rules.size(); idx++)
	{
		a[idx] = possible_reaction_rules[idx].k();
		const ReactionRule::reactant_container_type &reactants = 
			possible_reaction_rules[idx].reactants();
		for( ReactionRule::reactant_container_type::iterator it = reactants.begin();
				it != reactants.end();
				it++ )
		{
			a[idx] *= this->world_->num_molecules(*it);
		}
	}
	double a_total( std::accumulate(a.begin(), a.end(), double(0.0)) );
	if(a_total == 0.0)
	{
		return;
	}

	double rnd_num1(this->rng_.uniform(0, 1));
	double dt( gsl_sf_log(1.0 / rnd_num1) / double(a_total) );
	double rnd_num2(this->rng_.uniform(0, 1) * a_total);

	int u(-1);
	double acc(0.0);
	int len( a.size() );
	do
	{
		u++;
		acc += a[u];
	} while(acc < rnd_num2 && u < len -1);
	
	//Reaction[u] occurs.
	for( 
		ReactionRule::reactant_container_type::iterator it(possible_reaction_rules[u].reactants().begin());
		it != possible_reaction_rules[u].reactants().end();
		it++ )
	{
		int one(1);
		this->world_->remove_molecules(*it, one);
	}
	for(
		ReactionRule::product_container_type::iterator it(possible_reaction_rules[u].products().begin());
		it != possible_reaction_rules[u].products().end();
		it++ )
	{
		int one(1);
		this->world_->add_molecules(*it, one);
	}
	this->world_->set_t( this->world_->t() + dt );
	this->num_steps_++;
	//GillespieSolver gs(*(this->model_), *(this->world_), this->rng_);
	//gs.step();
}

bool GillespieSimulator::step(Real const &upto) 
{
	while(this->world_->t() < upto)
	{
		this->step();
	}
	return true;
}

void GillespieSimulator::run(void) 
{
	;
}

void GillespieSimulator::set_t(Real const &t) 
{
	this->world_->set_t(t);
}

Real GillespieSimulator::t(void) const 
{
	return this->world_->t();
}

Integer GillespieSimulator::num_steps(void) const
{
	return this->num_steps_;
}

RandomNumberGenerator& GillespieSimulator::rng(void) 
{
	return this->rng_;
}


}	// gillespie

}	// ecell4

