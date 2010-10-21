#ifndef PARTICLE_SIMULATION_STRUCTURE_HPP
#define PARTICLE_SIMULATION_STRUCTURE_HPP

#include "Structure.hpp"

template<typename Ttraits_>
struct ImmutativeStructureVisitor;

template<typename Ttraits_>
struct MutativeStructureVisitor;

template<typename Ttraits_>
struct ParticleSimulationStructure: public Structure<typename Ttraits_::world_type::traits_type>
{
    typedef Ttraits_ traits_type;

    virtual ~ParticleSimulationStructure() {}

    virtual void accept(ImmutativeStructureVisitor<traits_type> const&) const = 0;

    virtual void accept(MutativeStructureVisitor<traits_type> const&) = 0;
};

#endif /* PARTICLE_SIMULATION_STRUCTURE_HPP */
