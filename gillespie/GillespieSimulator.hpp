#ifndef __GILLESPIESIMULATOR_HPP
#define __GILLESPIESIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <core/RandomNumberGenerator.hpp>
#include <core/Model.hpp>
#include <core/NetworkModel.hpp>
#include <core/Simulator.hpp>

#include "GillespieWorld.hpp"
//#include "GillespieSolver.hpp"

namespace ecell4 
{

namespace gillespie 
{

class GillespieSimulator 
	: 
		public Simulator 
{
public:
	GillespieSimulator(
		boost::shared_ptr<NetworkModel> model, boost::shared_ptr<GillespieWorld> world,
		RandomNumberGenerator &rng)
		: model_(model), world_(world), rng_(rng)
	{
		this->num_steps_ = 0;
	}
		
	Integer num_steps(void) const;
	void step(void) ;
	bool step(Real const & upto);
	void run(void);

	Real t(void) const;
	void set_t(Real const &t);

	RandomNumberGenerator &rng(void);
					
protected:
	boost::shared_ptr<NetworkModel> model_;
	boost::shared_ptr<GillespieWorld> world_;
	
	Integer num_steps_;
	RandomNumberGenerator &rng_;
};

}

}	// ecell4

#endif //__GILLESPIESIMULATOR_HPP
