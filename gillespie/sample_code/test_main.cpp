#include <string>
#include <sstream>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "../GillespieSimulator.hpp"

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

	////////////////////////////////////////
	const int RANK = 1;
	hsize_t dim[] = {2};
	DataSpace space(RANK, dim);
	////////////////////////////////////////

	// ここで一覧をかく
	H5File *file = new H5File("CompartmentSpace.h5", H5F_ACC_TRUNC);
	std::string species_list_dataset = "species";
	typedef struct specie_name {
		uint32_t id;
		char name[32];
	} specie_name;
	specie_name list[2];
	CompType comp_specie( sizeof(specie_name) );
	std::string string_id("id");
	std::string string_name("name");
	comp_specie.insertMember( string_id, HOFFSET(specie_name, id), PredType::STD_I32LE);
	comp_specie.insertMember( string_name, HOFFSET(specie_name, name), StrType(PredType::C_S1, 32) );

	list[0].id = 1;
	strcpy(list[0].name, sp1.name().c_str());
	list[1].id = 2;
	strcpy(list[1].name, sp2.name().c_str());
	DataSpace space0(RANK, dim);

	DataSet *dataset_species = new DataSet(file->createDataSet(std::string("species"), comp_specie, space0));
	dataset_species->write(list, comp_specie);

	//グループは時間ごとのデータの集合
	// 構造体一つあたりの情報を登録してあげる
	typedef struct species_set {
		uint32_t id;
		uint32_t num;
	} species_set;
	species_set ls[2];
	CompType mtype1( sizeof(species_set) );
	std::string ids = "species_id";
	std::string numbers = "number";
	mtype1.insertMember( ids, HOFFSET(species_set, id), PredType::STD_I32LE);
	mtype1.insertMember( numbers, HOFFSET(species_set, num), PredType::STD_I32LE);

	for(int i = 0; i < 10; i++) {
		sim.step();

		std::ostringstream ost;
		ost << sim.t();
		DataSet *dataset = new DataSet(file->createDataSet(std::string( ost.str() ), mtype1, space));
		ls[0].id = 1;
		ls[0].num = world->num_molecules(sp1);
		ls[1].id = 2;
		ls[1].num = world->num_molecules(sp2);

		dataset->write(ls, mtype1);
		printf("t = %f A: %llu, B: %llu \n", sim.t(), world->num_molecules(sp1), world->num_molecules(sp2));
		ost.clear();

		delete dataset;
	}
	delete file;

}
