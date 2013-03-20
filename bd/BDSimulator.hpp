#ifndef __ECELL4_BD_BD_SIMULATOR_HPP
#define __ECELL4_BD_BD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/Simulator.hpp>

#include "BDWorld.hpp"
#include "BDPropagator.hpp"

#include <H5Cpp.h>
#include <hdf5.h>


namespace ecell4
{

namespace bd
{

class BDSimulator
    : public Simulator
{
public:

    BDSimulator(boost::shared_ptr<Model> model, boost::shared_ptr<BDWorld> world)
        : model_(model), world_(world), dt_(0), num_steps_(0)
    {
        ;
    }

    // SimulatorTraits

    Real t() const
    {
        return (*world_).t();
    }

    Real dt() const
    {
        return dt_;
    }

    Integer num_steps() const
    {
        return num_steps_;
    }

    void step();
    bool step(const Real& upto);

    // Optional members

    void set_t(const Real& t)
    {
        (*world_).set_t(t);
    }

    void set_dt(const Real& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

	void save_hdf5_init(std::string);
	void save_hdf5(void);

protected:

    boost::shared_ptr<Model> model_;
    boost::shared_ptr<BDWorld> world_;


	H5::H5File *file_;
	typedef struct h5_particles {
		int h5_particle_id;
		double h5_particle_position[3];
	} h5_particles;

	typedef struct h5_particles_index {
		int h5_particle_id;
		char h5_particle_name[32];
	} h5_particles_index;

    /**
     * the protected internal state of BDSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */
    Real dt_;
    Integer num_steps_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_SIMULATOR_HPP */
