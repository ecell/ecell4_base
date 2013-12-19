#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include "LatticeWorld.hpp"

namespace ecell4
{

namespace lattice
{

class LatticeSimulator
    : public Simulator
{
public:

    LatticeSimulator(
            boost::shared_ptr<NetworkModel> model,
            boost::shared_ptr<LatticeWorld> world)
        : model_(model), world_(world)
    {
        // TODO
    }

    virtual Real t() const
    {
        return (*world_).t();
    }

    virtual Real dt() const
    {
        return dt_;
    }

    Integer num_steps() const
    {
        return (Integer)(t() / dt());
    }

    void initialize();
    void step();
    bool step(const Real& upto);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<LatticeWorld> world_;

    Real dt_;
};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
