#ifndef __PARTICLE_SPACE_HPP
#define __PARTICLE_SPACE_HPP

#include <cmath>
#include <gsl/gsl_pow_int.h>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "functions.hpp"
#include "exceptions.hpp"
#include "Position3.hpp"
#include "Particle.hpp"
#include "Species.hpp"
#include "Space.hpp"

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

const H5std_string DATASET_NAME( "TimePoint" );
const H5std_string MEMBER1( "particle_id" );
const H5std_string MEMBER2( "positions" );

namespace ecell4
{

class ParticleSpace
    : public Space
{
public:

    typedef std::vector<std::pair<ParticleID, Particle> >
    particle_container_type;

public:

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
        throw NotImplemented("num_particles(const Species&) not implemented");
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
        throw NotImplemented("list_particles(const Species&) not implemented");
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

    // ParticleSpace member functions

    /**
     * update a particle specified by its ID.
     * if the particle does not exist, create a new particle.
     * this function is a member of ParticleSpace
     * @param pid ParticleID
     * @param p Particle
     * @return if the particle already exists or not bool
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
                retval += gsl_pow_2(diff - edge_length);
            }
            else if (diff < -half)
            {
                retval += gsl_pow_2(diff + edge_length);
            }
            else
            {
                retval += gsl_pow_2(diff);
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

    virtual void save(const std::string& filename) = 0;
};

class ParticleSpaceVectorImpl
    : public ParticleSpace
{
public:

    typedef ParticleSpace::particle_container_type particle_container_type;

protected:

    typedef utils::get_mapper_mf<
    ParticleID, particle_container_type::size_type>::type particle_map_type;

public:

    ParticleSpaceVectorImpl(const Position3& edge_lengths)
    {
        set_edge_lengths(edge_lengths);
    }

    // ParticleSpaceTraits

    const Position3& edge_lengths() const
    {
        return edge_lengths_;
    }

    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const;
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

    // Optional members

    const particle_container_type& particles() const
    {
        return particles_;
    }

    void save(const std::string& filename)
    {
        typedef struct h5_partcles {
            int h5_particle_id;
            double h5_particle_position[3];
        } h5_particles;

        h5_particles h5_p[particles_.size()];

        for (unsigned int i(0); i < particles_.size(); ++i)
        {
            h5_p[i].h5_particle_id = particles_[i].first;
            h5_p[i].h5_particle_position[0] = particles_[i].second.position()[0];
            h5_p[i].h5_particle_position[1] = particles_[i].second.position()[1];
            h5_p[i].h5_particle_position[2] = particles_[i].second.position()[2];
        }

        H5::Exception::dontPrint();

        // const H5std_string FILE_NAME("hoge.h5");
        H5File* file = new H5File(filename, H5F_ACC_RDONLY);
        CompType mtype(sizeof(h5_particles));
        mtype.insertMember(MEMBER1, HOFFSET(h5_particles, h5_particle_id),
                           PredType::NATIVE_INT);
        mtype.insertMember(MEMBER2, HOFFSET(h5_particles, h5_particle_position),
                           PredType::NATIVE_DOUBLE);

        hsize_t dim[] = {2, particles_.size()};
        DataSpace space(1, dim);

        DataSet* dataset;
        dataset = new DataSet(file->createDataSet(DATASET_NAME, mtype, space));
        dataset->write(h5_p, mtype);

        delete dataset;
        delete file;
    }

private:

    void set_edge_lengths(const Position3& edge_lengths);

protected:

    Position3 edge_lengths_;
    particle_container_type particles_;
    particle_map_type index_map_;
};

} // ecell4

#endif /* __PARTICLE_SPACE_HPP */
