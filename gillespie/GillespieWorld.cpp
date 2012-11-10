#include <map>
#include <vector>
#include <string>
#include <pficommon/text/json.h>

#include "./GillespieWorld.hpp"

using namespace pfi::text::json;

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

World *init_world_from_json(pfi::text::json::json js) {
	
	return new World();
}


#ifdef unit_world
int main(void)
{
	using namespace std;
	json json_world;
	ifstream ifs("./init.json");
	string f_str;
	
	ifs.seekg(0, ios::end);
	int len(ifs.tellg());
	ifs.seekg(0, ios::beg);
	
	char *buffer = new char[len + 1];
	ifs.read(buffer, len);
	stringstream ss(string(buffer));

	delete[] buffer;

}

#endif
