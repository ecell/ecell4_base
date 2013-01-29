#include "GillespieSimulator.hpp"
#include <numeric>
#include <vector>
#include <gsl/gsl_sf_log.h>
#include <sstream>

namespace ecell4 
{

namespace gillespie 
{

void GillespieSimulator::calc_next_reaction_(void) 
{
	// reset
	this->dt_ = inf;

	const NetworkModel::reaction_rule_container_type &possible_reaction_rules =
		this->model_->reaction_rules();
	if (possible_reaction_rules.size() == 0)
	{
		this->dt_ = inf;
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
		// Any reactions cannot occur.
		this->dt_ = inf;
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

	if (len == u)
	{
		// Any reactions cannot occur.
		this->dt_ = inf;
		return;
	}

	// save.
	this->next_reaction_num_ = u;
	this->dt_ = dt;
}

void GillespieSimulator::step(void)
{
	if (this->dt_ == inf) 
	{
		// Any reactions cannot occur.
		return;
	}
	const NetworkModel::reaction_rule_container_type &possible_reaction_rules =
		this->model_->reaction_rules();

	int u = this->next_reaction_num_;
	Real dt = this->dt_;

	if (dt == 0.0 || u < 0)
	{
		// Any reactions cannot occur.
		return;
	}

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

	this->calc_next_reaction_();	
}

bool GillespieSimulator::step(Real const &upto) 
{
	// proceed reactions before the argument 'upto'.
	while(this->dt_ != inf && this->world_->t() + this->dt_ < upto) 
	{
		this->step();
	}

	// The next reaction will occur after the argument 'upto'.
	if (this->dt_ != inf)
	{
		this->dt_ = this->t() + this->dt_ - upto ;
	}
	this->set_t(upto);
	return true;
}

void GillespieSimulator::initialize(void)
{
	this->calc_next_reaction_();
}

void GillespieSimulator::set_t(Real const &t) 
{
	this->world_->set_t(t);
}

Real GillespieSimulator::t(void) const 
{
	return this->world_->t();
}

Real GillespieSimulator::dt(void) const
{
	return this->dt_;
}

Integer GillespieSimulator::num_steps(void) const
{
	return this->num_steps_;
}

RandomNumberGenerator& GillespieSimulator::rng(void) 
{
	return this->rng_;
}

void GillespieSimulator::save_hdf5_init(std::string filename)
{
	using namespace H5;
	this->file_ = new H5File(filename, H5F_ACC_TRUNC);

	// save species' id
	typedef struct specie_id_struct {
		uint32_t id;
		char name[32];
	} specie_id_struct;

	const NetworkModel::species_container_type &species_list = this->model_->species();

	CompType mtype_specie_id( sizeof(specie_id_struct) );

	mtype_specie_id.insertMember(std::string("id"), HOFFSET(specie_id_struct, id), 
					PredType::STD_I32LE);
	mtype_specie_id.insertMember(std::string("name"), HOFFSET(specie_id_struct, name), 
					StrType(PredType::C_S1, 32));

	specie_id_struct *specie_id_table = new specie_id_struct[ species_list.size() ];
	for(unsigned int i = 0; i < species_list.size(); i++) 
	{
		specie_id_table[i].id = i + 1;
		strcpy(specie_id_table[i].name, species_list[i].name().c_str());
	}
	const int RANK = 1;
	hsize_t dim[1];	// id, name
	dim[0] = species_list.size();
	DataSpace space(RANK, dim);
	DataSet *dataset_species = new DataSet(this->file_->createDataSet(std::string("species"),
							mtype_specie_id, space));
	dataset_species->write(specie_id_table, mtype_specie_id);

	delete specie_id_table;
}

void GillespieSimulator::save_hdf5(void) 
{
	using namespace H5;
	typedef struct species_num_struct {
		uint32_t id;
		uint32_t num;
	} species_num_struct;

	// Construct Datatype.
	CompType mtype_species_num(sizeof(species_num_struct));
	mtype_species_num.insertMember(std::string("id"), HOFFSET(species_num_struct, id),
					PredType::STD_I32LE);
	mtype_species_num.insertMember(std::string("number"), HOFFSET(species_num_struct, num),
					PredType::STD_I32LE);

	const NetworkModel::species_container_type &species_list = this->model_->species();
	species_num_struct *species_num_table = new species_num_struct[ species_list.size() ];

	// Construct Data Set.
	for(unsigned int i = 0; i < species_list.size(); i++)
	{
		species_num_table[i].id = i + 1;
		species_num_table[i].num = this->world_->num_molecules( species_list[i] );
	}
	const int RANK = 1;
	hsize_t dim[1];
	dim[0] = species_list.size();
	std::ostringstream ost;
	ost << this->t();
	DataSpace space(RANK, dim);
	DataSet *dataset = new DataSet(this->file_->createDataSet(std::string( ost.str() ), 
							mtype_species_num, space));
	dataset->write(species_num_table, mtype_species_num);
	ost.clear();
	
	delete species_num_table;
}

}	// gillespie

}	// ecell4

