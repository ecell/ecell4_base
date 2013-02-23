#ifndef __ECELL4_SPACE_HPP
#define __ECELL4_SPACE_HPP

#include <stdexcept>

#include "exceptions.hpp"
#include "types.hpp"
#include "Position3.hpp"
#include "Species.hpp"
#include "Particle.hpp"


#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

namespace ecell4
{

class Space
{
public:

    Space()
        : t_(0.0)
    {
        ;
    }

    // SpaceTraits

    const Real& t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    // CompartmentSpaceTraits

    /**
     * get volume.
     * this function is a part of the trait of CompartmentSpace.
     * @return a volume (m^3) Real
     */
    virtual const Real& volume() const
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

    // ParticleSpaceTraits

    /**
     * get the axes lengths of a cuboidal region.
     * this function is a part of the trait of ParticleSpace.
     * @return edge lengths Position3
     */
    virtual const Position3& edge_lengths() const
    {
        throw NotSupported(
            "edge_lengths() is not supported by this space class");
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

protected:

    Real t_;
};

} // ecell4

#endif /* __ECELL4_SPACE_HPP */
