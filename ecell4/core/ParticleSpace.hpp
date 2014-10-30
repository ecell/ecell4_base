#ifndef __ECELL4_PARTICLE_SPACE_HPP
#define __ECELL4_PARTICLE_SPACE_HPP

#include <cmath>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "functions.hpp"
#include "exceptions.hpp"
#include "Position3.hpp"
#include "Particle.hpp"
#include "Species.hpp"
#include "Space.hpp"
#include "ParticleSpaceHDF5Writer.hpp"


namespace ecell4
{

class ParticleSpace
    : public Space
{
public:

    typedef std::vector<std::pair<ParticleID, Particle> >
    particle_container_type;

public:

    ParticleSpace()
        : t_(0.0)
    {
        ;
    }

    virtual ~ParticleSpace()
    {
        ; // do nothing
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

    // ParticleSpaceTraits

    /**
     * get the axes lengths of a cuboidal region.
     * this function is a part of the trait of ParticleSpace.
     * @return edge lengths Position3
     */
    virtual const Position3& edge_lengths() const
    {
        throw NotImplemented("edge_lengths() not implemented");
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @return a number of particles Integer
     */
    virtual Integer num_particles() const
    {
        throw NotImplemented("num_particles() not implemented");
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @param sp a species
     * @return a number of particles Integer
     */
    virtual Integer num_particles(const Species& sp) const
    {
        return num_particles_exact(sp);
    }

    virtual Integer num_particles_exact(const Species& sp) const
    {
        throw NotImplemented("num_particles_exact(const Species&) not implemented");
    }

    /**
     * get all particles.
     * this function is a part of the trait of ParticleSpace.
     * @return a list of particles
     */
    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles() const
    {
        throw NotImplemented("list_particles() not implemented");
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
        return list_particles_exact(sp);
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const
    {
        throw NotImplemented("list_particles_exact(const Species&) not implemented");
    }

    /**
     * check if the particle exists.
     * this function is a part of the trait of ParticleSpace.
     * @param pid an ID for the particle
     * @return if the particle exists or not bool
     */
    virtual bool has_particle(const ParticleID& pid) const
    {
        throw NotImplemented("has_particle(const ParticleID&) not implemented.");
    }

    virtual void save(H5::Group* root) const = 0;
    virtual void load(const H5::Group& root) = 0;

    // ParticleSpace member functions

    /**
     * update a particle specified by its ID.
     * if the particle does not exist, create a new particle.
     * this function is a member of ParticleSpace
     * @param pid ParticleID
     * @param p Particle
     * @return if the particle does not exist or not bool
     */
    virtual bool update_particle(const ParticleID& pid, const Particle& p) = 0;

    /**
     * get a particle specified with an ID.
     * this function is a member of ParticleSpace
     * @param pid ParticleID
     * @return a pair of ParticleID and Particle
     */
    virtual std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const = 0;

    /**
     * remove a particle
     * this function is a member of ParticleSpace
     * @param pid ParticleID
     */
    virtual void remove_particle(const ParticleID& pid) = 0;

    /**
     * get particles within a spherical region.
     * this function is a part of the trait of ParticleSpace.
     * @param pos a center position of the sphere
     * @param radius a radius of the sphere
     * @return a list of particles
     */
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius) const = 0;

    /**
     * get particles within a spherical region except for ignore(s).
     * this function is a part of the trait of ParticleSpace.
     * @param pos a center position of the sphere
     * @param radius a radius of the sphere
     * @param ignore an ignored ID
     * @return a list of particles
     */
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius,
        const ParticleID& ignore) const = 0;

    /**
     * get particles within a spherical region except for ignore(s).
     * this function is a part of the trait of ParticleSpace.
     * @param pos a center position of the sphere
     * @param radius a radius of the sphere
     * @param ignore1 an ignored ID
     * @param ignore2 an ignored ID
     * @return a list of particles
     */
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const = 0;

    /**
     * transpose a position based on the periodic boundary condition.
     * this function is a part of the trait of ParticleSpace.
     * @param pos1 a target position
     * @param pos2 a reference position
     * @return a transposed position Position3
     */
    Position3 periodic_transpose(
        const Position3& pos1, const Position3& pos2) const
    {
        Position3 retval(pos1);
        const Position3& edges(edge_lengths());
        for (Position3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edges[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

            if (diff > half)
            {
                retval[dim] += edge_length;
            }
            else if (diff < -half)
            {
                retval[dim] -= edge_length;
            }
        }
        return retval;
    }

    /**
     * transpose a position based on the periodic boundary condition.
     * if the position is in the region, returns the original position.
     * this function is a part of the trait of ParticleSpace.
     * @param pos a target position
     * @return a transposed position Position3
     */
    inline Position3 apply_boundary(const Position3& pos) const
    {
        return modulo(pos, edge_lengths());
    }

    /**
     * calculate a square of the distance between positions
     * this function is a part of the trait of ParticleSpace.
     * @param pos1
     * @param pos2
     * @return a square of the distance
     */
    Real distance_sq(
        const Position3& pos1, const Position3& pos2) const
    {
        Real retval(0);
        const Position3& edges(edge_lengths());
        for (Position3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edges[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

            if (diff > half)
            {
                retval += pow_2(diff - edge_length);
            }
            else if (diff < -half)
            {
                retval += pow_2(diff + edge_length);
            }
            else
            {
                retval += pow_2(diff);
            }
        }
        return retval;
    }

    /**
     * calculate the distance between positions
     * this function is a part of the trait of ParticleSpace.
     * @param pos1
     * @param pos2
     * @return the distance
     */
    inline Real distance(const Position3& pos1, const Position3& pos2) const
    {
        return std::sqrt(distance_sq(pos1, pos2));
    }

    // Optional members

    virtual const particle_container_type& particles() const = 0;

    virtual Real get_value(const Species& sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

protected:

    Real t_;
};

class ParticleSpaceVectorImpl
    : public ParticleSpace
{
public:

    typedef ParticleSpace base_type;
    typedef ParticleSpace::particle_container_type particle_container_type;

protected:

    typedef utils::get_mapper_mf<
    ParticleID, particle_container_type::size_type>::type particle_map_type;

public:

    ParticleSpaceVectorImpl(const Position3& edge_lengths)
    {
        reset(edge_lengths);
    }

    // ParticleSpaceTraits

    const Position3& edge_lengths() const
    {
        return edge_lengths_;
    }

    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const;
    bool has_particle(const ParticleID& pid) const;

    // ParticleSpace member functions

    bool update_particle(const ParticleID& pid, const Particle& p);
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const;
    void remove_particle(const ParticleID& pid);

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius,
        const ParticleID& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const;

    // CompartmentSpaceTraits

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    // Optional members

    const particle_container_type& particles() const
    {
        return particles_;
    }

    void save(H5::Group* root) const
    {
        save_particle_space(*this, root);
    }

    void load(const H5::Group& root)
    {
        load_particle_space(root, this);
    }

    void reset(const Position3& edge_lengths);

protected:

    Position3 edge_lengths_;
    particle_container_type particles_;
    particle_map_type index_map_;
};

} // ecell4

#endif /* __ECELL4_PARTICLE_SPACE_HPP */
