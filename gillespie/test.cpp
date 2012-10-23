#include "./Gillespie.hpp"

int main(void)
{
	GillespieSolver gs;
	int world[] = {1000, 1000, 1000, 1000};
	gs.set_current_state(world, sizeof(world)/sizeof(int));

#define TEMP_ID(x)	x-'W'
	int ri1 = gs.reaction_add();
	gs.reaction_add_substance(ri1, TEMP_ID('X'), 1);
	gs.reaction_add_product(ri1, TEMP_ID('Y'), 1);
	gs.reaction_set_kinetic_parameter(ri1, 0.5);

	int ri2 = gs.reaction_add();
	gs.reaction_add_substance(ri2, TEMP_ID('Y'), 1);
	gs.reaction_add_product(ri2, TEMP_ID('X'), 1);
	gs.reaction_set_kinetic_parameter(ri2, 0.2);

	int ri3 = gs.reaction_add();
	gs.reaction_add_substance(ri3, TEMP_ID('X'), 2);
	gs.reaction_add_product(ri3, TEMP_ID('Z'), 1);
	gs.reaction_set_kinetic_parameter(ri3, 0.4);

	int ri4 = gs.reaction_add();
	gs.reaction_add_substance(ri4, TEMP_ID('Z'), 1);
	gs.reaction_add_product(ri4, TEMP_ID('X'), 2);
	gs.reaction_set_kinetic_parameter(ri4, 0.2);

	int ri5 = gs.reaction_add();
	gs.reaction_add_substance(ri5, TEMP_ID('X'), 1);
	gs.reaction_add_substance(ri5, TEMP_ID('W'), 1);
	gs.reaction_add_product(ri5, TEMP_ID('X'), 2);
	gs.reaction_set_kinetic_parameter(ri5, 0.3);

	int ri6 = gs.reaction_add();
	gs.reaction_add_substance(ri6, TEMP_ID('X'), 2);
	gs.reaction_add_product(ri6, TEMP_ID('X'), 1);
	gs.reaction_add_product(ri6, TEMP_ID('W'), 1);
	gs.reaction_set_kinetic_parameter(ri6, 0.5);

	for(double elapse = 0.0; elapse < 10.0; elapse += gs.run(1.0)) {
		fprintf(stderr,
		"%f, %d, %d, %d, %d\n",
		gs.get_current_time(),
		gs.current_state[1],
		gs.current_state[2],
		gs.current_state[3],
		gs.current_state[0]	);
	}
	return 0;
}
