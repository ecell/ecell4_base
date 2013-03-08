#ifndef __BD_SIMULATOR_HPP
#define __BD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Simulator.hpp>

#include "BDWorld.hpp"
#include "BDPropagator.hpp"

#include <hdf5.h>
#include <H5Cpp.h>


namespace ecell4
{

namespace bd
{

class BDSimulator
    : public Simulator
{
public:

    BDSimulator(
        boost::shared_ptr<Model> model, boost::shared_ptr<BDWorld> world,
        RandomNumberGenerator& rng)
        : model_(model), world_(world), rng_(rng), num_steps_(0), dt_(0)
    {
    	// about hdf5
    	this->file_ = NULL;
        ;
    }
    ~BDSimulator(void)
    {
    	if (this->file_ != NULL)
    	{
    		delete this->file_;
    	}
    }


    Real t() const
    {
        return (*world_).t();
    }

    void set_t(Real const& t)
    {
        (*world_).set_t(t);
    }

    // about hdf5
    void save_hdf5_init(std::string filename);
    void save(void);

    Real dt() const
    {
        return dt_;
    }

    void set_dt(Real const& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

    Integer num_steps() const
    {
        return num_steps_;
    }

    RandomNumberGenerator& rng()
    {
        return rng_;
    }

    void step();
    bool step(Real const& upto);


protected:

    boost::shared_ptr<Model> model_;
    boost::shared_ptr<BDWorld> world_;

    // about hdf5
    H5::H5File *file_;

    /**
     * the protected internal state of BDSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */
    Real dt_;
    Integer num_steps_;
    RandomNumberGenerator& rng_;
};

} // bd

} // ecell4

#endif /* __BD_SIMULATOR_HPP */
