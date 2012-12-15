#include "GillespieSimulator.hpp"

namespace ecell4 
{

namespace gillespie 
{

void GillespieSimulator::step(void) 
{
	GillespieSolver gs(*(this->model_), *(this->world_), this->rng_);
	gs.step();
}

bool GillespieSimulator::step(Real const &upto) 
{
	printf("%s\n", __func__);
	return true;
}

void GillespieSimulator::run(void) 
{
	;
}

void GillespieSimulator::set_t(Real const &t) 
{
	;
}

Real GillespieSimulator::t(void) const 
{
	return this->world_->t();
}

Integer GillespieSimulator::num_steps(void) const
{
	return this->num_steps_;
}

RandomNumberGenerator& GillespieSimulator::rng(void) 
{
	return this->rng_;
}


}	// gillespie

}	// ecell4

