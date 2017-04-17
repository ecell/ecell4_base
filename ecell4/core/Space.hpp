#ifndef ECELL4_SPACE_HPP
#define ECELL4_SPACE_HPP

#include <stdexcept>

#include "exceptions.hpp"
#include "types.hpp"
#include "Real3.hpp"
#include "Species.hpp"
#include "Particle.hpp"

#include <iostream>

namespace ecell4
{

class Space
{
public:

    typedef enum
    {
        ELSE,
        PARTICLE,
        LATTICE,
        COMPARTMENT,
        SUBVOLUME
    } space_kind;

public:

    Space()
    {
        ;
    }

    virtual ~Space()
    {
        ; // do nothing
    }

    // SpaceTraits

    // virtual const Real t() const = 0;  //XXX: This doesn't work with Python35 MSVC2015
    virtual const Real t() const
    {
        return 0.0;  //XXX: Just for debugging
    }

    virtual void set_t(const Real& t) = 0;

    virtual void save(const std::string& filename) const = 0;
    //XXX: This doesn't work with Python35 MSVC2015
    // virtual void save(const std::string& filename) const
    // {
    //     throw NotSupported(
    //         "save(const std::string) is not supported by this space class");
    // }

    virtual void load(const std::string& filename)
    {
        throw NotSupported(
            "load(const std::string) is not supported by this space class");
    }

    // CompartmentSpaceTraits

    /**
     * get volume.
     * this function is a part of the trait of CompartmentSpace.
     * @return a volume (m^3) Real
     */
    virtual const Real volume() const
    {
        throw NotSupported("volume() is not supported by this space class");
    }

    /**
     * get the number of species in this space.
     * this function is a part of the trait of CompartmentSpace.
     * @return a number of species Integer
     */
    virtual Integer num_species() const
    {
        throw NotSupported("num_species() is not supported by this space class");
    }

    /**
     * return if the species is in this space or not.
     * this function is a part of the trait of CompartmentSpace.
     * @param sp a species
     * @return if the species is in this space
     */
    virtual bool has_species(const Species& sp) const
    {
        throw NotSupported(
            "has_species(const Species&) is not supported by this space class");
    }

    /**
     * get the number of molecules
     * this function is a part of the trait of CompartmentSpace.
     * @param sp a species
     * @return a number of molecules Integer
     */
    virtual Integer num_molecules(const Species& sp) const
    {
        throw NotSupported(
            "num_molecules(const Species&) is not supported"
            " by this space class");
    }

    virtual Integer num_molecules_exact(const Species& sp) const
    {
        throw NotSupported(
            "num_molecules_exact(const Species&) is not supported"
            " by this space class");
    }

    virtual Real get_value(const Species& sp) const
    {
        throw NotSupported(
            "get_value(const Species&) is not supported"
            " by this space class");
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        throw NotSupported(
            "get_value_exact(const Species&) is not supported"
            " by this space class");
    }

    // ParticleSpaceTraits

    /**
     * get the axes lengths of a cuboidal region.
     * this function is a part of the trait of ParticleSpace.
     * @return edge lengths Real3
     */
    virtual const Real3& edge_lengths() const
    {
        throw NotSupported(
            "edge_lengths() is not supported by this space class");
    }

    /**
     * get the actual axes lengths of a cuboidal region.
     * this function is a part of the trait of ParticleSpace.
     * return edge lengths as a default.
     * overload this function if the actual size is not equal to
     * edge lengths.
     * @return actual edge lengths Real3
     */
    virtual Real3 actual_lengths() const
    {
        return edge_lengths();
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @return a number of particles Integer
     */
    virtual Integer num_particles() const
    {
        throw NotSupported(
            "num_particles() is not supported by this space class");
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @param sp a species
     * @return a number of particles Integer
     */
    virtual Integer num_particles(const Species& sp) const
    {
        throw NotSupported(
            "num_particles(const Species&) is not supported"
            " by this space class");
    }

    virtual Integer num_particles_exact(const Species& sp) const
    {
        throw NotSupported(
            "num_particles_exact(const Species&) is not supported"
            " by this space class");
    }

    /**
     * check if the particle exists.
     * this function is a part of the trait of ParticleSpace.
     * @param pid an ID for the particle
     * @return if the particle exists or not bool
     */
    virtual bool has_particle(const ParticleID& pid) const
    {
        throw NotSupported(
            "has_particle(const ParticleID&) is not supported"
            " by this space class");
    }

    virtual std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {
        throw NotSupported(
            "get_particle(const ParticleID&) is not supported"
            " by this space class");
    }

    /**
     * get all particles.
     * this function is a part of the trait of ParticleSpace.
     * @return a list of particles
     */
    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles() const
    {
        throw NotSupported(
            "list_particles() is not supported by this space class.");
    }

    /**
     * get particles.
     * this function is a part of the trait of ParticleSpace.
     * @param sp a species
     * @return a list of particles
     */
    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const
    {
        throw NotSupported(
            "list_particles(const Species&) is not supported"
            " by this space class");
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const
    {
        throw NotSupported(
            "list_particles_exact(const Species&) is not supported"
            " by this space class");
    }
};

} // ecell4

#endif /* ECELL4_SPACE_HPP */
