#include "./GillespieWorld.hpp"
#include "./GillespieSolver.hpp"

int main(void)
{
	World my_world;
	Model my_model;

#define TEMP_ID(x)	x-'W'
	ReactionRule r1;
	r1.add_reactant(TEMP_ID('X'), 1);
	r1.add_product(TEMP_ID('Y'), 1);
	r1.set_kinetic_parameter(0.5);
	my_model.reactions.push_back(r1);

	ReactionRule r2;
	r2.add_reactant(TEMP_ID('Y'), 1);
	r2.add_product(TEMP_ID('X'), 1);
	r2.set_kinetic_parameter(0.2);
	my_model.reactions.push_back(r2);

	ReactionRule r3;
	r3.add_reactant(TEMP_ID('X'), 2);
	r3.add_product(TEMP_ID('Z'), 1);
	r3.set_kinetic_parameter(0.4);
	my_model.reactions.push_back(r3);

	ReactionRule r4;
	r4.add_reactant(TEMP_ID('Z'), 1);
	r4.add_product(TEMP_ID('X'), 2);
	r4.set_kinetic_parameter(0.2);
	my_model.reactions.push_back(r4);
	
	ReactionRule r5;
	r5.add_reactant(TEMP_ID('X'), 1);
	r5.add_reactant(TEMP_ID('W'), 1);
	r5.add_product(TEMP_ID('X'), 2);
	r5.set_kinetic_parameter(0.3);
	my_model.reactions.push_back(r5);

	ReactionRule r6;
	r6.add_reactant(TEMP_ID('X'), 2);
	r6.add_product(TEMP_ID('X'), 1);
	r6.add_product(TEMP_ID('W'), 1);
	r6.set_kinetic_parameter(0.5);
	my_model.reactions.push_back(r6);

	my_world.add_specie(TEMP_ID('X'), 1000);
	my_world.add_specie(TEMP_ID('Y'), 1000);
	my_world.add_specie(TEMP_ID('Z'), 1000);
	my_world.add_specie(TEMP_ID('W'), 1000);

	GillespieSolver gs(my_world, my_model);
	for(double elapse = 0.0; elapse < 10.0; elapse += gs.run(1.0)) {
		fprintf(
				stderr,
				"%f, %d, %d, %d, %d\n",
				my_world.get_current_time(),
				my_world.get_current_state(TEMP_ID('X') ),
				my_world.get_current_state(TEMP_ID('Y') ),
				my_world.get_current_state(TEMP_ID('Z') ),
				my_world.get_current_state(TEMP_ID('W') )	);
	}
	return 0;
}

