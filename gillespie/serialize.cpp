
#include "pficommon/text/json.h"
#include "pficommon/text/csv.h"
#include "GillespieWorld.hpp"
#include "GillespieSolver.hpp"
#include <cstdlib>


using namespace pfi::text::json;
using namespace std;
//============================================================
//	Common
//============================================================

json string_to_json(std::string str) {
	json js;
	std::stringstream ss(str);
	ss >> js;
	return js;
}

string read_file_all(const char *json_filename) {
	ifstream ifs(json_filename);
	string content;
	if (ifs) {
		ifs.seekg(0, ios::end);
		int len(ifs.tellg());
		ifs.seekg(0, ios::beg);
		char *buffer = new char[len];

		ifs.read(buffer, len);
		content = string(buffer);
		delete[] buffer;
	}
	return content;
}

//============================================================
//	World
//============================================================
World *init_world_from_json(json js_world) {
	World *world = new World();
	for(unsigned int idx = 0; idx < js_world.size(); idx++) {
		string species(json_cast<string>(js_world[idx]["species"]));
		int initVal(json_cast<int>(js_world[idx]["initVal"]));

		world->add_specie(species, initVal);
	}
	return world;
}


World *init_world_from_csv(string csv_str) {
	bool header = true, first = true;
	World *world = new World();
	pfi::text::csv_parser psr(csv_str);
	for(pfi::text::csv_iterator p(psr), q; p != q; ++p) {
		if (header == true && first == true) {
			first = false;
			continue;
		}
		string species( (*p)[0] );
		int initVal( atoi((*p)[1]) );
		world->add_specie(species, initVal);
	}
	return world;
}


//============================================================
//	Model
//============================================================
ReactionRule *init_reaction_from_json(json js_reaction) {
	ReactionRule *r = new ReactionRule;
	string reactant( json_cast<string>(js_reaction["reactant"]) );
	int reaStoich( json_cast<int>(js_reaction["reaStoich"]) );
	r->add_reactant(reactant, reaStoich);

	string product( json_cast<string>(js_reaction["product"]) );
	int proStoich( json_cast<int>(js_reaction["proStoich"]) );
	r->add_product(product, proStoich);
	r->set_kinetic_parameter( json_cast<double>(js_reaction["kineticParameter"]) );

	return r;
}


Model *init_model_from_json(json js_model) {
	Model *m = new Model;
	for(unsigned int i = 0; i < js_model.size(); i++) {
		ReactionRule *r = init_reaction_from_json(js_model[i]);
		m->reactions.push_back(*r);
		delete r;
	}
	return m;
}

Model *init_model_from_csv(string &csv_model) {
	bool header = true, first = true;
	Model *m = new Model;
	pfi::text::csv_parser psr(csv_model);
	for(pfi::text::csv_iterator p(psr), q; p != q; ++p) {
		if (header == true && first == true) {
			first = false;
			continue;
		}
		ReactionRule *r = new ReactionRule;
		string reactant( (*p)[1] );
		string product( (*p)[3] );
		r->add_reactant( reactant, atoi((*p)[2]) );
		r->add_product(product, atoi((*p)[4]) );
		r->set_kinetic_parameter( double(atof( (*p)[5] )) );
		m->reactions.push_back( *r );
		delete r;
	}
	return m;
}


