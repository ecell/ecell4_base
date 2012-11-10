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
