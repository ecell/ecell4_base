
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
template<typename F>
World *init_world_from_json(json js_world, F translate_func ) {
	World *world = new World();
	for(unsigned int idx = 0; idx < js_world.size(); idx++) {
		string species(json_cast<string>(js_world[idx]["species"]));
		int initVal(json_cast<int>(js_world[idx]["initVal"]));

		world->add_specie(translate_func(species), initVal);
	}
	return world;
}


template <typename F>
World *init_world_from_csv(string csv_str, F translate_func) {
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
		world->add_specie(translate_func(species), initVal);
	}
	return world;
}


//============================================================
//	Model
//============================================================
template <typename F>
ReactionRule *init_reaction_from_json(json js_reaction, F translate_func) {
	ReactionRule *r = new ReactionRule;
	string reactant( json_cast<string>(js_reaction["reactant"]) );
	int reaStoich( json_cast<int>(js_reaction["reaStoich"]) );
	r->add_reactant(translate_func(reactant), reaStoich);

	string product( json_cast<string>(js_reaction["product"]) );
	int proStoich( json_cast<int>(js_reaction["proStoich"]) );
	r->add_product(translate_func(product), proStoich);
	r->set_kinetic_parameter( json_cast<double>(js_reaction["kineticParameter"]) );

	return r;
}

template <typename F>
Model *init_model_from_json(json js_model, F translate_func) {
	Model *m = new Model;
	for(unsigned int i = 0; i < js_model.size(); i++) {
		ReactionRule *r = init_reaction_from_json(js_model[i], translate_func);
		m->reactions.push_back(*r);
		delete r;
	}
	return m;
}

//Specieと内部におけるそのIDとの変換をどうしようか。。。？
//とりあえずはこのファンクタを変換テーブルの代わりにした。
class Specie_to_Id {
public:
	int operator()(std::string &specie) {
		int id;
		if (specie == "X") {
			id = 1;
		} else if (specie == "Y") {
			id = 2;
		} else if (specie == "Z") {
			id = 3;
		} else {
			id = 4;
		}
		return id;
	}
};

#ifdef unit_world
int main(void)
{
	Specie_to_Id translater;
	//string json_file_content(read_file_all("./data/init.json"));
	//World *w = init_world_from_json( string_to_json(json_file_content), translater );
	string csv_file_content(read_file_all("./data/init.csv"));
	World *w = init_world_from_csv( csv_file_content, translater );
	std::cout << *w << std::endl;

	string json_model(read_file_all("./data/reaction.json"));
	Model *m = init_model_from_json( string_to_json(json_model), translater );

	GillespieSolver gs(*w, *m);
	for(double elapse = 0.0; elapse < 10.0; elapse += gs.run(1.0)) {
		std::cout << *w << std::endl;
	}
}
#endif
