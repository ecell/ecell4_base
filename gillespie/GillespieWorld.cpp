#include <map>
#include <vector>
#include <string>
#include <pficommon/text/json.h>

#include "./GillespieWorld.hpp"


//============================================================
//	World	*Definitions
//============================================================
World::World(void)
{
	this->current_t = 0.0;
}
double World::get_current_time(void)
{	return this->current_t;	}

void World::set_current_time(double new_t)
{	this->current_t = new_t;	}

int World::get_current_state(int id)
{	
	return this->current_state[id];	
}

void World::set_current_state(int id, int number)
{	
	this->current_state[id] = number;	
}

void World::add_specie(int id, int number = 0)
{
	this->current_state.insert(std::map<int,int>::value_type(id, number));
}

using namespace pfi::text::json;
using namespace std;
int spiecie_to_id(string specie)
{
	int id;
	if(specie == "X") {
		id = 1;
	} else if (specie == "Y") {
		id = 2;
	} else {
		id = 3;
	}
	return id;
}

// world -- factorial method
World *init_world_from_json(json js_world) {
	World *world = new World();
	for(int idx = 0; idx < js_world.size(); idx++) {
		string species(json_cast<string>(js_world[idx]["species"]));
		int initVal(json_cast<int>(js_world[idx]["initVal"]));
		world->add_specie(spiecie_to_id(species), initVal);
	}
	return world;
}

json string_to_json(std::string str) {
	json js;
	std::stringstream ss(str);
	ss >> js;
	return js;
}

// XXX should be core api ?
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

#ifdef unit_world
int main(void)
{
	string json_file_content(read_file_all("./data/init.json"));
	World *w = init_world_from_json( string_to_json(json_file_content) );

}
#endif
