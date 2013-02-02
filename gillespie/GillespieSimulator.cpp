#include "GillespieSimulator.hpp"
#include <numeric>
#include <vector>
#include <gsl/gsl_sf_log.h>
#include <sstream>

#include <boost/scoped_array.hpp>

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
}

void GillespieSimulator::save_hdf5(void)
{
	using namespace H5;
	if (this->file_ == NULL)
		return;

	// Define Data Structure
	CompType mtype_id_table_struct(sizeof(species_id_table_struct));
	mtype_id_table_struct.insertMember(std::string("id"), HOFFSET(species_id_table_struct, id),
					PredType::STD_I32LE);
	mtype_id_table_struct.insertMember(std::string("name"), 
					HOFFSET(species_id_table_struct, name), StrType(PredType::C_S1, 32));

	CompType mtype_num_struct(sizeof(species_num_struct));
	mtype_num_struct.insertMember(std::string("id"), HOFFSET(species_num_struct, id),
							PredType::STD_I32LE);
	mtype_num_struct.insertMember(std::string("number"), HOFFSET(species_num_struct, num_of_molecules),
							PredType::STD_I32LE);

	// Construct Data Set.
	const NetworkModel::species_container_type &species_list = this->model_->species();
	boost::scoped_array<species_id_table_struct> species_id_table(new species_id_table_struct[ species_list.size() ]);
	boost::scoped_array<species_num_struct> species_num_table(new species_num_struct[ species_list.size() ]);

	for(unsigned int i(0); i < species_list.size(); i++) 
	{
		species_id_table[i].id = i + 1;
		strcpy(species_id_table[i].name, species_list[i].name().c_str());

		species_num_table[i].id = i + 1;
		species_num_table[i].num_of_molecules = this->world_->num_molecules( species_list[i] );
	}
	const int RANK = 1;
	hsize_t dim[1];
	dim[0] = species_list.size();

	// Create Path.
	std::ostringstream ost_hdf5path;
	ost_hdf5path << "/" << this->t();
	boost::scoped_ptr<Group> group (new Group(this->file_->createGroup( ost_hdf5path.str() )));

	DataSpace space(RANK, dim);
	std::string species_table_path = ost_hdf5path.str() + "/species";
	std::string species_num_path = ost_hdf5path.str() + "/num";
	boost::scoped_ptr<DataSet> dataset_id_table( new DataSet(this->file_->createDataSet(species_table_path , mtype_id_table_struct, space)) );
	boost::scoped_ptr<DataSet> dataset_num_table( new DataSet(this->file_->createDataSet(species_num_path , mtype_num_struct, space)) );
	dataset_id_table->write(species_id_table.get(), mtype_id_table_struct);
	dataset_num_table->write(species_num_table.get(), mtype_num_struct);
}

}	// gillespie

}	// ecell4

