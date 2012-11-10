#include <map>
#include <vector>

//============================================================
//	World	
//============================================================
class World {
private:
public:	// Data members
	double current_t;
	std::map<int,int> current_state;	// [id] -> number of substrate	// id => 分子の個数　

public:
	World();
	double get_current_time(void);
	void set_current_time(double new_t);
	int get_current_state(int id);
	void set_current_state(int id, int number);
	void add_specie(int id, int number);
};
