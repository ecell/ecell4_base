#include <pficommon/text/json.h>

#ifndef INCLUDE_GUARD_GILLESPIE_WORLD
#	define INCLUDE_GUARD_GILLESPIE_WORLD
//============================================================
//	World	
//============================================================
class World {
private:
public:	// Data members
	double current_t;
	std::map<std::string,int> current_state;	// [id] -> number of substrate	// id => 分子の個数　

public:
	World();
	double get_current_time(void);
	void set_current_time(double new_t);
	int get_current_state(std::string &sp);
	void set_current_state(std::string &sp, int number);
	void add_specie(std::string &sp, int number);

	std::string to_string(void);
};

template<typename F>
World *init_world_from_json(pfi::text::json::json js, F translate_func);

std::ostream &operator<<(std::ostream &s, World &w);

#endif	//INCLUDE_GUARD_GILLESPIE_WORLD
