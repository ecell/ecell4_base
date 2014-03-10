#ifndef PARTICLE_SIMULATOR_FACTORY_HPP
#define PARTICLE_SIMULATOR_FACTORY_HPP

#include <boost/noncopyable.hpp>
#include "ParticleSimulator.hpp"
#include "ParticleModel.hpp"

template<typename Ttraits_>
class ParticleSimulatorFactory
{
public:
    typedef Ttraits_ traits_type;

public:
    ParticleSimulatorFactory() {}

    virtual ~ParticleSimulatorFactory() {}

    virtual ParticleSimulator<traits_type>* operator()(ParticleModel const& model) const = 0;
};

#endif /* PARTICLE_SIMULATION_HPP */
