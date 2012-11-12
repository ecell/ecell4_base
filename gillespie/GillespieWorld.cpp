#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
#include "./GillespieWorld.hpp"

//============================================================
//	World	
//		-- General
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

//============================================================
//	World
//		-- pretty printer
//============================================================
string World::to_string(void) {
	ostringstream os;
	os << "time: " << this->current_t;
	for(std::map<int,int>::iterator it = this->current_state.begin(); it != this->current_state.end(); it++) {
		os << ", " << it->first << ": " << it->second;
	}
	return os.str();
}

ostream &operator<<(ostream &s, World &w) {
	return s << w.to_string();
}

