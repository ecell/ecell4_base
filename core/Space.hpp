#ifndef __SPACE_HPP
#define __SPACE_HPP

#include <stdexcept>

#include "exceptions.hpp"
#include "types.hpp"
#include "Position3.hpp"
#include "Species.hpp"
#include "Particle.hpp"


namespace ecell4
{

class Space
{
public:

    Space()
        : t_(0)
    {
        ;
    }

    Real const& t() const
    {
        return t_;
    }

    void set_t(Real const& t)
    {
        if (t < 0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    /**
     * get volume.
     * this function is a part of the trait of CompartmentSpace.
     * @return a volume (m^3) Real
     */
    virtual Real const& volume() const
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
    virtual bool has_species(Species const& sp) const
    {
        throw NotSupported("has_species() is not supported by this space class");
    }

    /**
     * get the number of molecules
     * this function is a part of the trait of CompartmentSpace.
     * @param sp a species
     * @return a number of molecules Integer
     */
    virtual Integer num_molecules(Species const& sp) const
    {
        throw NotSupported("num_molecules() is not supported by this space class");
    }

    /**
     * get the axes lengths of a cuboidal region.
     * this function is a part of the trait of ParticleSpace.
     * @return edge lengths Position3
     */
    virtual Position3 const& edge_lengths() const
    {
        throw NotSupported("edge_lengths() is not supported by this space class");
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @return a number of particles Integer
     */
    virtual Integer num_particles() const
    {
        throw NotSupported("num_particles() is not supported by this space class");
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @param sp a species
     * @return a number of particles Integer
     */
    virtual Integer num_particles(Species const& species) const
    {
        throw NotSupported("num_particles() is not supported by this space class");
    }

    /**
     * get particles.
     * this function is a part of the trait of ParticleSpace.
     * @param sp a species
     * @return a list of particles
     */
    virtual std::vector<std::pair<ParticleID, Particle> >
    get_particles(Species const& species) const
    {
        throw NotSupported("get_particles() is not supported by this space class");
    }

protected:

    Real t_;
};

} // ecell4

#endif /* __SPACE_HPP */
