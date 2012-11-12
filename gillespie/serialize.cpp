
#include "pficommon/text/json.h"
#include "GillespieWorld.hpp"
#include "GillespieSolver.hpp"


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


//============================================================
//	Model
//============================================================
ReactionRule *init_reaction_from_json(json js_reaction) {
	ReactionRule *r = new ReactionRule;
	return r;
}

Model *init_model_from_json(json js_model) {
	Model *m = new Model;
	for(unsigned int i = 0; i < js_model.size(); i++) {
		m->reactions.push_back( *init_reaction_from_json(js_model[i]) );
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
	string json_file_content(read_file_all("./data/init.json"));
	World *w = init_world_from_json( string_to_json(json_file_content), translater );
	std::cout << *w << std::endl;

	string json_model(read_file_all("./data/reactions.json"));

}
#endif
