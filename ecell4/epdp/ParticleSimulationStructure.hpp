#ifndef PARTICLE_SIMULATION_STRUCTURE_HPP
#define PARTICLE_SIMULATION_STRUCTURE_HPP

#include "Structure.hpp"

template<typename Ttraits_>
struct ImmutativeStructureVisitor;

template<typename Ttraits_>
struct MutativeStructureVisitor;

// template<typename Ttraits_>
// struct ParticleSimulationStructure: public Structure<typename Ttraits_::world_type::traits_type>
// {
//     typedef Ttraits_ traits_type;
//     typedef Structure<typename traits_type::world_type::traits_type> base_type;
//     typedef typename base_type::identifier_type identifier_type;
// 
//     virtual ~ParticleSimulationStructure() {}
// 
//     virtual void accept(ImmutativeStructureVisitor<traits_type> const&) const = 0;
// 
//     virtual void accept(MutativeStructureVisitor<traits_type> const&) = 0;
// 
//     ParticleSimulationStructure(identifier_type const& id): base_type(id) {}
// };

template<typename Ttraits_>
struct ParticleSimulationStructure: public Structure<Ttraits_>
{
    typedef Ttraits_ traits_type;
    typedef Structure<traits_type> base_type;
    typedef typename base_type::identifier_type identifier_type;

    virtual ~ParticleSimulationStructure() {}

    virtual void accept(ImmutativeStructureVisitor<traits_type> const&) const = 0;

    virtual void accept(MutativeStructureVisitor<traits_type> const&) = 0;

    ParticleSimulationStructure(identifier_type const& id): base_type(id) {}
};

#endif /* PARTICLE_SIMULATION_STRUCTURE_HPP */
